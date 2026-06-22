#!/usr/bin/env python3
"""
ga_optimize_halbach.py — otimizador genetico (Cooley et al. 2018) para o Halbach
tubular esparso do CS-MPFM-LFMRI. Fecha o campo em B0_alvo e minimiza a
nao-homogeneidade (ppm) no FOV do tubo, respeitando bore minimo e massa maxima.

Genes (vetor de projeto p):
  r1, r2        raios das 2 camadas concentricas [m]
  n1, n2        contagem de imas por camada no CORPO
  n1e, n2e      contagem por camada nas SLICES de EXTREMIDADE (end-compensation)
  grade         indice {0:N52 BR=1.48, 1:N42 BR=1.32}

Topologia fixa: 9 slices (0,+-20,+-40,+-60,+-80 mm), 2 extremidades end-comp/lado.
Fitness vetorizado: reproduz EXATAMENTE a matematica dipolar de HalbachRing.calculateB
(B em Tesla, momento BR*dim^3, fator 1/4pi, rotacao k=2). Sem gmsh/solid (numpy puro).

Saida: JSON no schema do HalbachMRIDesigner, pronto p/ '--fem'.
"""
import argparse
import json
import numpy as np

GRADES = [("N52", 1.48), ("N42", 1.32)]
DIM = 0.015          # cubo 15 mm
MUR = 1.05
RHO_NDFEB = 7500.0   # kg/m^3
SLICE_SEP = 0.020
N_SLICES = 9         # posicoes -80..80 mm
END_IDX = {0, 1, N_SLICES - 2, N_SLICES - 1}  # 2 extremidades por lado


def ring_B(radius, N, dim, BR, zpos, ev, k=2):
    """Campo dipolar de 1 anel (replica de HalbachRing.calculateB), B em Tesla."""
    ang = np.linspace(0, 2 * np.pi, N, endpoint=False)
    x, y, z = ev
    B = np.zeros((x.size, 3))
    for a in ang:
        px, py = radius * np.cos(a), radius * np.sin(a)
        mabs = BR * dim ** 3
        ma = a * k
        m2 = np.array([mabs * np.cos(ma), mabs * np.sin(ma)]) / (4 * np.pi)
        rx, ry, rz = x + px, y + py, z + zpos
        rdm = 3 * (rx * m2[0] + ry * m2[1])
        r2 = rx * rx + ry * ry + rz * rz
        r3 = r2 ** 1.5
        r5 = r2 ** 2.5
        B[:, 0] += rx * rdm / r5 - m2[0] / r3
        B[:, 1] += ry * rdm / r5 - m2[1] / r3
        B[:, 2] += rz * rdm / r5
    return B


PACK = 1.35   # folga azimutal: cubo girado (k=2) tem pegada ate ~dim*sqrt(2)
RADGAP = 1.10 # folga radial minima entre camadas (centro-a-centro >= dim*RADGAP)


def n_max(radius):
    """Maximo de cubos 'dim' que cabem no perimetro sem colidir (passo angular >= dim*PACK)."""
    return int(np.floor(2 * np.pi * radius / (DIM * PACK)))


def evaluate(gene, ev, center_mask):
    r1, r2, n1, n2, n1e, n2e, grade = gene
    BR = GRADES[int(round(grade))][1]
    # restricao de empacotamento: cubos nao podem colidir no anel
    nm1, nm2 = n_max(r1), n_max(r2)
    overflow = (max(0, int(n1) - nm1) + max(0, int(n2) - nm2)
                + max(0, int(n1e) - nm1) + max(0, int(n2e) - nm2))
    n1, n1e = min(int(n1), nm1), min(int(n1e), nm1)
    n2, n2e = min(int(n2), nm2), min(int(n2e), nm2)
    zs = (np.arange(N_SLICES) - (N_SLICES - 1) / 2) * SLICE_SEP
    B = np.zeros((ev[0].size, 3))
    n_imas = 0
    for i, z0 in enumerate(zs):
        a, b = (n1e, n2e) if i in END_IDX else (n1, n2)
        B += ring_B(r1, a, DIM, BR, z0, ev)
        B += ring_B(r2, b, DIM, BR, z0, ev)
        n_imas += a + b
    mag = np.linalg.norm(B[:, :2], axis=1)
    Bc = mag[center_mask].mean()
    ppm = (mag.max() - mag.min()) / mag.mean() * 1e6
    mass = n_imas * DIM ** 3 * RHO_NDFEB
    bore = 2 * (r1 - DIM / 2)
    if (r2 - r1) < DIM * RADGAP:                 # colisao radial entre camadas
        overflow += 100
    return Bc, ppm, mass, bore, n_imas, overflow


def fitness(gene, ev, center_mask, target_B0, bore_min, mass_max):
    Bc, ppm, mass, bore, _, overflow = evaluate(gene, ev, center_mask)
    J = abs(Bc - target_B0) / target_B0          # erro de campo (dominante)
    J += 0.02 * ppm / 1e4                          # homogeneidade
    if bore < bore_min:
        J += 50 * (bore_min - bore) / bore_min     # restricao bore
    if mass > mass_max:
        J += 5 * (mass - mass_max) / mass_max      # restricao massa
    J += 0.5 * overflow                            # restricao empacotamento (colisao)
    return J


def make_roi(r_fov=0.020, z_fov=0.025, nr=9, nz=11):
    xs = np.linspace(-r_fov, r_fov, nr)
    zs = np.linspace(-z_fov, z_fov, nz)
    pts = [(x, y, z) for x in xs for y in xs for z in zs if x * x + y * y <= r_fov * r_fov]
    pts = np.array(pts)
    ev = [pts[:, 0].copy() + 1e-6, pts[:, 1].copy(), pts[:, 2].copy()]
    cmask = (np.abs(pts[:, 0]) < 1e-3) & (np.abs(pts[:, 1]) < 1e-3) & (np.abs(pts[:, 2]) < 1e-3)
    if not cmask.any():
        cmask = np.linalg.norm(pts, axis=1) < r_fov / 3
    return ev, cmask


# limites dos genes: r1, r2, n1, n2, n1e, n2e, grade
LO = np.array([0.048, 0.068, 12, 18, 14, 20, 0])
HI = np.array([0.064, 0.090, 24, 34, 28, 38, 1.4])


def run_ga(target_B0, bore_min, mass_max, pop, gen, seed=0):
    rng = np.random.default_rng(seed)
    ev, cmask = make_roi()
    P = LO + (HI - LO) * rng.random((pop, LO.size))
    best, best_J = None, np.inf
    for g in range(gen):
        J = np.array([fitness(ind, ev, cmask, target_B0, bore_min, mass_max) for ind in P])
        order = np.argsort(J)
        P = P[order]
        J = J[order]
        if J[0] < best_J:
            best_J, best = J[0], P[0].copy()
        if g % 10 == 0 or g == gen - 1:
            Bc, ppm, mass, bore, ni, _ = evaluate(P[0], ev, cmask)
            print(f"  gen {g:3d}  J={J[0]:.4f}  B0={Bc*1000:6.1f}mT  ppm={ppm:6.0f}  "
                  f"massa={mass:4.1f}kg  bore={bore*1000:4.1f}mm  N={ni}")
        # elitismo + crossover + mutacao
        elite = P[: max(2, pop // 5)]
        children = []
        while len(children) < pop - elite.shape[0]:
            pa, pb = elite[rng.integers(len(elite))], elite[rng.integers(len(elite))]
            mask = rng.random(LO.size) < 0.5
            child = np.where(mask, pa, pb)
            child += rng.normal(0, 0.06, LO.size) * (HI - LO)  # mutacao gaussiana
            child = np.clip(child, LO, HI)
            children.append(child)
        P = np.vstack([elite, np.array(children)])
    return best, ev, cmask


def to_json(gene):
    r1, r2, n1, n2, n1e, n2e, grade = gene
    gi = int(round(grade))
    name, BR = GRADES[gi]
    nm1, nm2 = n_max(r1), n_max(r2)
    n1, n1e = min(int(round(n1)), nm1), min(int(round(n1e)), nm1)
    n2, n2e = min(int(round(n2)), nm2), min(int(round(n2e)), nm2)
    rings = [
        {"id": 0, "radius": round(r1 * 1e3, 1), "numMagnets": n1, "_role": "interna-corpo"},
        {"id": 1, "radius": round(r2 * 1e3, 1), "numMagnets": n2, "_role": "externa-corpo"},
        {"id": 2, "radius": round(r1 * 1e3, 1), "numMagnets": n1e, "_role": "interna-extremidade"},
        {"id": 3, "radius": round(r2 * 1e3, 1), "numMagnets": n2e, "_role": "externa-extremidade"},
    ]
    positions = [0, 20, 40, 60, 80]
    slices = []
    for p in positions:
        end = p >= 60
        sl = {"position": str(p), "innerRadius": "45", "outerRadius": "95",
              "rings": [{"id": 2 if end else 0}, {"id": 3 if end else 1}]}
        if p not in (0, 80):
            sl["standWidth"] = "0"
            sl["standHeight"] = "0"
        slices.append(sl)
    return {
        "_comment": f"GA-otimizado (Cooley 2018). Grade {name} BR={BR}. Cubo {DIM*1e3:.0f}mm.",
        "magnets": [{"_grade": name, "dimension": str(DIM * 1e3), "shape": "cube",
                     "BR": str(BR), "mur": str(MUR)}],
        "defaultMagnetType": 0,
        "rings": rings,
        "mirrorSlices": True,
        "shimTrayHeight": "14", "shimTrayAngle": "25", "shimTrayRadius": "90",
        "standWidth": "200", "standHeight": "150",
        "numConnectionRods": 8, "connectionRodsArcRadius": 100, "connectionRodsDiameter": 5,
        "slices": slices,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("seed_json", nargs="?", default="examples/mpfm_03T_dn50.json")
    ap.add_argument("--target-B0", type=float, default=0.30)
    ap.add_argument("--bore-min", type=float, default=0.090)
    ap.add_argument("--mass-max", type=float, default=12.0)
    ap.add_argument("--pop", type=int, default=40)
    ap.add_argument("--gen", type=int, default=80)
    ap.add_argument("-o", default="examples/mpfm_03T_opt.json")
    args = ap.parse_args()

    print(f"[GA] alvo B0={args.target_B0} T  bore>={args.bore_min*1e3:.0f}mm  "
          f"massa<={args.mass_max}kg  pop={args.pop} gen={args.gen}")
    best, ev, cmask = run_ga(args.target_B0, args.bore_min, args.mass_max, args.pop, args.gen)
    Bc, ppm, mass, bore, ni, _ = evaluate(best, ev, cmask)
    print(f"[GA] MELHOR: B0={Bc*1000:.1f}mT (f0={42.577*Bc*1e3/1e3:.3f}MHz)  ppm={ppm:.0f}  "
          f"massa={mass:.1f}kg  bore={bore*1e3:.1f}mm  N={ni}  grade={GRADES[int(round(best[6]))][0]}")
    cfg = to_json(best)
    with open(args.o, "w") as f:
        json.dump(cfg, f, indent=4, ensure_ascii=False)
    print(f"[GA] escrito: {args.o}")
