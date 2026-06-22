#!/usr/bin/env python3
"""
ga_homog.py — GA com objetivo de homogeneidade SHIMMAVEL (versao 2 CAMADAS).
Esta e' a versao que produziu o BOM resultado: unshim=599 / raw=900 / B0=294,6 mT
(seed 0, pop 60, gen 150, b0-lo 0.27). Bandas de raio APERTADAS por camada (a
generalizacao p/ N camadas afrouxou isso e degradou a convergencia -> ver
ga_homog_multilayer.py, mantido a parte).

Objetivo (chave): minimizar o RESIDUO INALCANCAVEL pelo shim (p-p de Δppm - proj_{SH<=3})
SUJEITO a o conteudo shimmavel caber na capacidade real do shim (raw <= SHIM_CAP),
pois o shim tem limite de corrente (nao adianta dumpar amplitude impossivel no
subespaco shimmavel). Topologia: 5 grupos-z (0,20,40,60,80mm) espelhados -> 9 slices,
2 camadas. Genes (14): r1,r2 | n[z][layer] (5x2=10) | g1,g2 (grade por camada).
"""
import argparse, json
import numpy as np

GRADES = [("N52", 1.48), ("N42", 1.32)]
DIM = 0.015; MUR = 1.05; RHO = 7500.0
SEP = 0.020; ZPOS = np.array([0, 20, 40, 60, 80]) * 1e-3
PACK = 1.35; RADGAP = 1.10   # 1.35 = otimo campo+homog (294mT). PACK 1.40 testado: cai p/ ~270mT (caro).
# Colisao de mesh do escolhido = artefato OCC (cubos tocando) -> resolver encolhendo caixa no FEM, nao tirando ima.
ROI_R = 0.018; ROI_Z = 0.025
SHIM_CAP = 2500.0   # capacidade realista do shim (ppm p-p) — passivo+ativo c/ limite de corrente


def ring_B(radius, N, BR, zpos, ev, k=2):
    if N < 1: return 0.0
    ang = np.linspace(0, 2*np.pi, int(N), endpoint=False)
    x, y, z = ev; B = np.zeros((x.size, 3))
    for a in ang:
        px, py = radius*np.cos(a), radius*np.sin(a)
        mabs = BR*DIM**3; ma = a*k
        m2 = np.array([mabs*np.cos(ma), mabs*np.sin(ma)])/(4*np.pi)
        rx, ry, rz = x+px, y+py, z+zpos
        rdm = 3*(rx*m2[0]+ry*m2[1]); r2 = rx*rx+ry*ry+rz*rz
        r3 = r2**1.5; r5 = r2**2.5
        B[:, 0] += rx*rdm/r5 - m2[0]/r3
        B[:, 1] += ry*rdm/r5 - m2[1]/r3
        B[:, 2] += rz*rdm/r5
    return B


def n_max(r):
    return int(np.floor(2*np.pi*r/(DIM*PACK)))


def shim_basis(pts):
    """Harmonicos solidos cartesianos ate ordem 3 = subespaco SHIMMAVEL (passivo+ativo).
    Inclui ordem 0 (constante=offset/lock). Coluna por modo, normalizadas."""
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    cols = [np.ones_like(x),
            x, y, z,                                              # ordem 1
            2*z*z-x*x-y*y, x*z, y*z, x*x-y*y, x*y,                # ordem 2
            z*(2*z*z-3*x*x-3*y*y), x*(4*z*z-x*x-y*y), y*(4*z*z-x*x-y*y),
            z*(x*x-y*y), x*y*z, x*(x*x-3*y*y), y*(3*x*x-y*y)]     # ordem 3
    A = np.stack(cols, axis=1).astype(float)
    A /= np.linalg.norm(A, axis=0, keepdims=True)
    return A


def make_roi(nr=9, nz=9):
    xs = np.linspace(-ROI_R, ROI_R, nr); zs = np.linspace(-ROI_Z, ROI_Z, nz)
    pts = np.array([(x, y, z) for x in xs for y in xs for z in zs if x*x+y*y <= ROI_R*ROI_R])
    ev = [pts[:, 0].copy()+1e-6, pts[:, 1].copy(), pts[:, 2].copy()]
    cmask = (np.abs(pts[:, 0]) < 1e-3) & (np.abs(pts[:, 1]) < 1e-3) & (np.abs(pts[:, 2]) < 1e-3)
    if not cmask.any(): cmask = np.linalg.norm(pts, axis=1) < ROI_R/3
    Bsh = shim_basis(pts)
    Pproj = Bsh @ np.linalg.pinv(Bsh)        # projetor no subespaco shimmavel
    return ev, cmask, Pproj


def decode(gene):
    r1, r2 = gene[0], gene[1]
    counts = gene[2:12].reshape(5, 2)
    g1, g2 = int(round(gene[12])), int(round(gene[13]))
    return r1, r2, counts, g1, g2


def evaluate(gene, ev, cmask, Pproj):
    r1, r2, counts, g1, g2 = decode(gene)
    BR1, BR2 = GRADES[g1][1], GRADES[g2][1]
    nm1, nm2 = n_max(r1), n_max(r2)
    overflow = 0
    B = np.zeros((ev[0].size, 3)); n_imas = 0
    for zi, z0 in enumerate(ZPOS):
        a = min(int(round(counts[zi, 0])), nm1); b = min(int(round(counts[zi, 1])), nm2)
        a = max(a, 0); b = max(b, 0)
        overflow += max(0, int(round(counts[zi, 0]))-nm1) + max(0, int(round(counts[zi, 1]))-nm2)
        for zz in ({z0, -z0} if z0 > 0 else {0.0}):
            if a: B = B + ring_B(r1, a, BR1, zz, ev)
            if b: B = B + ring_B(r2, b, BR2, zz, ev)
            n_imas += a + b
    if (r2 - r1) < DIM*RADGAP: overflow += 100
    mag = np.linalg.norm(B[:, :2], axis=1)
    Bc = mag[cmask].mean()
    dppm = (mag - mag.mean())/mag.mean()*1e6
    ppm_raw = np.ptp(mag)/mag.mean()*1e6
    # RESIDUO INALCANCAVEL pelo shim = dppm menos projecao no subespaco shimmavel
    resid = dppm - Pproj @ dppm
    ppm_unshim = np.ptp(resid)
    mass = n_imas*DIM**3*RHO; bore = 2*(r1-DIM/2)
    return Bc, ppm_raw, ppm_unshim, mass, bore, n_imas, overflow


def fitness(gene, ev, cmask, Pproj, b0_lo, b0_hi, mass_max, bore_min):
    Bc, ppm_raw, ppm_unshim, mass, bore, ni, ov = evaluate(gene, ev, cmask, Pproj)
    # OBJETIVO: residuo INALCANCAVEL pelo shim, sem dumpar amplitude impossivel
    # no subespaco shimmavel (shim tem limite de corrente ~ SHIM_CAP ppm).
    J = ppm_unshim/1000.0
    if ppm_raw > SHIM_CAP: J += 5*(ppm_raw-SHIM_CAP)/SHIM_CAP   # raw alem do shimmavel
    if Bc < b0_lo: J += 30*(b0_lo-Bc)/b0_lo
    if Bc > b0_hi: J += 30*(Bc-b0_hi)/b0_hi
    if bore < bore_min: J += 50*(bore_min-bore)/bore_min
    if mass > mass_max: J += 5*(mass-mass_max)/mass_max
    J += 0.5*ov
    return J


# limites: r1,r2 | 10 counts | g1,g2. Bandas por-camada NAO-sobrepostas (convergencia),
# mas ALARGADAS p/ explorar mais campo. Recomputadas por DIM em set_bounds().
LO = np.array([0.050, 0.072] + [10]*10 + [0, 0])
HI = np.array([0.066, 0.094] + [26]*10 + [1.4, 1.4])

def set_bounds():
    global LO, HI
    r1_min = 0.045 + DIM/2 + 0.001          # bore >= 90mm
    # inner APERTADO no bore (campo alto); outer banda moderada. NAO alargar inner: o
    # objetivo de homog foge p/ raios grandes (campo fraco) se deixado solto.
    LO = np.array([r1_min,        r1_min+0.016] + [8]*10  + [0, 0])
    HI = np.array([r1_min+0.009,  r1_min+0.042] + [n_max(r1_min+0.042)]*10 + [1.4, 1.4])


def run(args):
    rng = np.random.default_rng(args.seed)
    ev, cmask, Pproj = make_roi()
    P = LO + (HI-LO)*rng.random((args.pop, LO.size))
    best, bestJ = None, np.inf
    n_imm = max(2, args.pop//10)                 # imigrantes aleatorios (anti-estagnacao)
    for g in range(args.gen):
        J = np.array([fitness(ind, ev, cmask, Pproj, args.b0_lo, args.b0_hi, args.mass_max, args.bore_min) for ind in P])
        order = np.argsort(J); P = P[order]; J = J[order]
        if J[0] < bestJ: bestJ, best = J[0], P[0].copy()
        if g % 15 == 0 or g == args.gen-1:
            Bc, raw, uns, mass, bore, ni, ov = evaluate(P[0], ev, cmask, Pproj)
            print(f"  gen {g:3d} J={J[0]:.3f} unshim={uns:5.0f} raw={raw:6.0f} B0={Bc*1e3:6.1f}mT "
                  f"massa={mass:4.1f}kg N={ni} ov={ov}")
        elite = P[:max(2, args.pop//5)]; kids = []
        while len(kids) < args.pop - elite.shape[0] - n_imm:
            pa, pb = elite[rng.integers(len(elite))], elite[rng.integers(len(elite))]
            child = np.where(rng.random(LO.size) < 0.5, pa, pb)
            child += rng.normal(0, 0.10, LO.size)*(HI-LO)
            kids.append(np.clip(child, LO, HI))
        imm = LO + (HI-LO)*rng.random((n_imm, LO.size))     # sangue novo
        P = np.vstack([elite, np.array(kids), imm])
    return best, ev, cmask, Pproj


def to_json(gene):
    r1, r2, counts, g1, g2 = decode(gene)
    nm1, nm2 = n_max(r1), n_max(r2)
    mags = [{"_grade": GRADES[g1][0], "dimension": str(DIM*1e3), "shape": "cube", "BR": str(GRADES[g1][1]), "mur": str(MUR)},
            {"_grade": GRADES[g2][0], "dimension": str(DIM*1e3), "shape": "cube", "BR": str(GRADES[g2][1]), "mur": str(MUR)}]
    rings, slices = [], []
    rid = 0
    for zi, z in enumerate(ZPOS*1e3):
        a = min(int(round(counts[zi, 0])), nm1); b = min(int(round(counts[zi, 1])), nm2)
        id_a, id_b = rid, rid+1; rid += 2
        rings.append({"id": id_a, "radius": round(r1*1e3, 1), "numMagnets": max(a, 1), "magnetType": 0})
        rings.append({"id": id_b, "radius": round(r2*1e3, 1), "numMagnets": max(b, 1), "magnetType": 1})
        sl = {"position": str(int(z)), "innerRadius": "45", "outerRadius": "95",
              "rings": [{"id": id_a}, {"id": id_b}]}
        if int(z) not in (0, 80): sl["standWidth"] = "0"; sl["standHeight"] = "0"
        slices.append(sl)
    return {"_comment": "GA-homogeneidade 2 camadas (objetivo: residuo shimmavel). grade-mista por camada.",
            "magnets": mags, "defaultMagnetType": 0, "rings": rings, "mirrorSlices": True,
            "shimTrayHeight": "14", "shimTrayAngle": "25", "shimTrayRadius": "90",
            "standWidth": "200", "standHeight": "150", "numConnectionRods": 8,
            "connectionRodsArcRadius": 100, "connectionRodsDiameter": 5, "slices": slices}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--b0-lo", type=float, default=0.27)
    ap.add_argument("--b0-hi", type=float, default=0.34)
    ap.add_argument("--bore-min", type=float, default=0.090)
    ap.add_argument("--mass-max", type=float, default=14.0)
    ap.add_argument("--pop", type=int, default=60)
    ap.add_argument("--gen", type=int, default=150)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--dim", type=float, default=0.015, help="aresta do cubo (m); maior -> mais campo (∝d²/anel)")
    ap.add_argument("-o", default="examples/mpfm_03T_homog.json")
    args = ap.parse_args()
    DIM = args.dim
    set_bounds()
    print(f"[GA-homog 2cam] objetivo=RESIDUO INALCANCAVEL pelo shim (SH<=3)  B0 in [{args.b0_lo},{args.b0_hi}]T  "
          f"ROI r<={ROI_R*1e3:.0f}mm  pop={args.pop} gen={args.gen} seed={args.seed}")
    best, ev, cmask, Pproj = run(args)
    Bc, raw, uns, mass, bore, ni, ov = evaluate(best, ev, cmask, Pproj)
    print(f"[GA-homog 2cam] MELHOR: unshim={uns:.0f}ppm  raw={raw:.0f}ppm  B0={Bc*1e3:.1f}mT (f0={42.577*Bc:.3f}MHz)  "
          f"massa={mass:.1f}kg bore={bore*1e3:.1f}mm N={ni}  grades=({GRADES[int(round(best[12]))][0]},{GRADES[int(round(best[13]))][0]})")
    json.dump(to_json(best), open(args.o, "w"), indent=2, ensure_ascii=False)
    print(f"[GA-homog 2cam] escrito: {args.o}")
