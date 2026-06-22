#!/usr/bin/env python3
"""
ga_homog.py — GA com objetivo de homogeneidade SHIMMAVEL, generalizado p/ L camadas
e tamanho de cubo configuravel (testar hipotese: cubos menores -> mais imas/anel ->
ordem azimutal mais alta -> menos p-p de borda; precisa de mais camadas p/ manter B0).

Objetivo (chave): minimizar o RESIDUO INALCANCAVEL pelo shim (p-p de Δppm - proj_{SH<=3})
SUJEITO a o conteudo shimmavel caber na capacidade real do shim (raw <= SHIM_CAP).
Genes: L raios + contagem por (grupo-z x camada) + grade por camada. Empacotamento
azimutal+radial respeitado. Anti-estagnacao com imigrantes.
"""
import argparse, json
import numpy as np

GRADES = [("N52", 1.48), ("N42", 1.32)]
DIM = 0.015; MUR = 1.05; RHO = 7500.0
ZPOS = np.array([0, 20, 40, 60, 80]) * 1e-3
NZ = len(ZPOS)
PACK = 1.35; RADGAP = 1.10
ROI_R = 0.018; ROI_Z = 0.025
SHIM_CAP = 2500.0
N_LAYERS = 2                       # sobrescrito por --layers
BORE = 0.092                       # bore livre alvo (m); sobrescrito por --bore (destrava o r1_min)
W_FIELD = 3.0                      # peso do premio de CAMPO ALTO (alvo 0,7T+, quanto maior melhor)


def ring_B(radius, N, BR, zpos, ev, k=2):
    # VETORIZADO (jun/2026): mesma fisica do loop por angulo, sem perda de qualidade.
    N = int(N)
    if N < 1: return 0.0
    ang = np.linspace(0, 2*np.pi, N, endpoint=False)          # (N,)
    x, y, z = ev                                              # (P,)
    px = radius*np.cos(ang); py = radius*np.sin(ang)          # (N,)
    mabs = BR*DIM**3; ma = ang*k
    m2x = mabs*np.cos(ma)/(4*np.pi); m2y = mabs*np.sin(ma)/(4*np.pi)   # (N,)
    rx = x[None, :] + px[:, None]                             # (N,P)
    ry = y[None, :] + py[:, None]
    rz = z[None, :] + zpos
    rdm = 3*(rx*m2x[:, None] + ry*m2y[:, None])
    r2 = rx*rx + ry*ry + rz*rz
    r3 = r2**1.5; r5 = r2**2.5
    Bx = (rx*rdm/r5 - m2x[:, None]/r3).sum(0)
    By = (ry*rdm/r5 - m2y[:, None]/r3).sum(0)
    Bz = (rz*rdm/r5).sum(0)
    return np.stack([Bx, By, Bz], axis=1)                     # (P,3)


def n_max(r): return int(np.floor(2*np.pi*r/(DIM*PACK)))


def shim_basis(pts):
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    cols = [np.ones_like(x), x, y, z,
            2*z*z-x*x-y*y, x*z, y*z, x*x-y*y, x*y,
            z*(2*z*z-3*x*x-3*y*y), x*(4*z*z-x*x-y*y), y*(4*z*z-x*x-y*y),
            z*(x*x-y*y), x*y*z, x*(x*x-3*y*y), y*(3*x*x-y*y)]
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
    return ev, cmask, Bsh @ np.linalg.pinv(Bsh)


def decode(gene):
    L = N_LAYERS
    radii = np.sort(gene[:L])
    counts = gene[L:L+NZ*L].reshape(NZ, L)
    grades = [int(round(g)) for g in gene[L+NZ*L:]]
    return radii, counts, grades


def evaluate(gene, ev, cmask, Pproj):
    radii, counts, grades = decode(gene)
    L = N_LAYERS
    nm = [n_max(r) for r in radii]
    overflow = 0
    B = np.zeros((ev[0].size, 3)); n_imas = 0
    for li in range(L):
        BR = GRADES[grades[li]][1]
        for zi, z0 in enumerate(ZPOS):
            a = int(round(counts[zi, li]))
            overflow += max(0, a - nm[li]); a = max(0, min(a, nm[li]))
            for zz in ({z0, -z0} if z0 > 0 else {0.0}):
                if a: B = B + ring_B(radii[li], a, BR, zz, ev)
                n_imas += a
    for li in range(1, L):                       # colisao radial entre camadas
        if radii[li]-radii[li-1] < DIM*RADGAP: overflow += 100
    mag = np.linalg.norm(B[:, :2], axis=1)
    Bc = mag[cmask].mean()
    dppm = (mag-mag.mean())/mag.mean()*1e6
    ppm_raw = np.ptp(mag)/mag.mean()*1e6
    ppm_unshim = np.ptp(dppm - Pproj @ dppm)
    mass = n_imas*DIM**3*RHO; bore = 2*(radii[0]-DIM/2)
    return Bc, ppm_raw, ppm_unshim, mass, bore, n_imas, overflow


def fitness(gene, ev, cmask, Pproj, b0_lo, b0_hi, mass_max, bore_min):
    Bc, raw, uns, mass, bore, ni, ov = evaluate(gene, ev, cmask, Pproj)
    J = uns/1000.0
    if raw > SHIM_CAP: J += 5*(raw-SHIM_CAP)/SHIM_CAP
    if Bc < b0_lo: J += 100*(b0_lo-Bc)/b0_lo   # piso -5% (b0_lo): penalidade FORTE
    J += W_FIELD*(b0_hi - min(Bc, b0_hi))/b0_hi # PREMIA campo alto (alvo 0,7T+); satura em b0_hi
    if bore < bore_min: J += 50*(bore_min-bore)/bore_min
    if mass > mass_max: J += 5*(mass-mass_max)/mass_max
    J += 0.5*ov
    return J


def bounds():
    # bandas de raio POR CAMADA, nao-sobrepostas, DIM-aware. inner APERTADO no bore (>=90mm);
    # passo radial entre camadas > DIM*RADGAP (sem colisao radial). banda de 10mm por camada.
    L = N_LAYERS
    r1_min = BORE/2 + DIM/2 + 0.001         # bore livre alvo = BORE (--bore); destravado do 92mm fixo
    step = DIM * 1.12                         # passo radial minimo (>RADGAP*DIM=22mm p/ 20mm) -> +campo/camada
    lo_r = [r1_min + i*step for i in range(L)]
    hi_r = [r1_min + 0.010 + i*step for i in range(L)]
    nmax = int(np.floor(2*np.pi*hi_r[-1]/(DIM*PACK)))
    lo = lo_r + [6]*(NZ*L) + [0]*L
    hi = hi_r + [nmax]*(NZ*L) + [1.4]*L
    return np.array(lo), np.array(hi)


def run(args):
    rng = np.random.default_rng(args.seed)
    ev, cmask, Pproj = make_roi()
    LO, HI = bounds()
    P = LO + (HI-LO)*rng.random((args.pop, LO.size))
    best, bestJ = None, np.inf
    n_imm = max(2, args.pop//10)
    stale = 0
    for g in range(args.gen):
        J = np.array([fitness(ind, ev, cmask, Pproj, args.b0_lo, args.b0_hi, args.mass_max, args.bore_min) for ind in P])
        order = np.argsort(J); P = P[order]; J = J[order]
        if J[0] < bestJ - 1e-9:
            bestJ, best = J[0], P[0].copy(); stale = 0
        else:
            stale += 1
        if stale >= args.patience:
            print(f"  early-stop gen {g} (plato {args.patience} gen sem melhora)"); break
        if g % 25 == 0 or g == args.gen-1:
            Bc, raw, uns, mass, bore, ni, ov = evaluate(P[0], ev, cmask, Pproj)
            print(f"  gen {g:3d} J={J[0]:.3f} unshim={uns:5.0f} raw={raw:6.0f} B0={Bc*1e3:6.1f}mT "
                  f"massa={mass:4.1f}kg N={ni} ov={ov}")
        elite = P[:max(2, args.pop//5)]; kids = []
        while len(kids) < args.pop - elite.shape[0] - n_imm:
            pa, pb = elite[rng.integers(len(elite))], elite[rng.integers(len(elite))]
            child = np.where(rng.random(LO.size) < 0.5, pa, pb)
            child += rng.normal(0, 0.10, LO.size)*(HI-LO)
            kids.append(np.clip(child, LO, HI))
        imm = LO + (HI-LO)*rng.random((n_imm, LO.size))
        P = np.vstack([elite, np.array(kids), imm])
    return best, ev, cmask, Pproj


def to_json(gene):
    radii, counts, grades = decode(gene)
    L = N_LAYERS; nm = [n_max(r) for r in radii]
    mags = [{"_grade": GRADES[grades[li]][0], "dimension": str(DIM*1e3), "shape": "cube",
             "BR": str(GRADES[grades[li]][1]), "mur": str(MUR)} for li in range(L)]
    rings, slices = [], []; rid = 0
    for zi, z in enumerate(ZPOS*1e3):
        sl_ring_ids = []
        for li in range(L):
            a = max(1, min(int(round(counts[zi, li])), nm[li]))
            rings.append({"id": rid, "radius": round(radii[li]*1e3, 1), "numMagnets": a, "magnetType": li})
            sl_ring_ids.append({"id": rid}); rid += 1
        sl = {"position": str(int(z)), "innerRadius": "44", "outerRadius": str(int(radii[-1]*1e3)+12),
              "rings": sl_ring_ids}
        if int(z) not in (0, 80): sl["standWidth"] = "0"; sl["standHeight"] = "0"
        slices.append(sl)
    return {"_comment": f"GA-homog L={L} camadas, cubo {DIM*1e3:.0f}mm. objetivo=residuo shimmavel.",
            "magnets": mags, "defaultMagnetType": 0, "rings": rings, "mirrorSlices": True,
            "shimTrayHeight": "14", "shimTrayAngle": "25", "shimTrayRadius": "90",
            "standWidth": "200", "standHeight": "150", "numConnectionRods": 8,
            "connectionRodsArcRadius": 100, "connectionRodsDiameter": 5, "slices": slices}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dim", type=float, default=0.015)
    ap.add_argument("--layers", type=int, default=2)
    ap.add_argument("--nz", type=int, default=5, help="nº de grupos-z (slices = 2*nz-1); axial. sep 20mm")
    ap.add_argument("--b0-lo", type=float, default=0.27)
    ap.add_argument("--b0-hi", type=float, default=0.34)
    ap.add_argument("--bore-min", type=float, default=0.088)
    ap.add_argument("--bore", type=float, default=0.0, help="bore livre alvo (m); destrava o r1_min (0=usa --bore-min/legado 92mm)")
    ap.add_argument("--mass-max", type=float, default=16.0)
    ap.add_argument("--pop", type=int, default=70)
    ap.add_argument("--gen", type=int, default=180)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--restarts", type=int, default=1, help="nº de sementes; guarda a melhor viável")
    ap.add_argument("--patience", type=int, default=80, help="early-stop: gerações sem melhora antes de parar")
    ap.add_argument("-o", default="examples/mpfm_03T_homog.json")
    args = ap.parse_args()
    DIM = args.dim; N_LAYERS = args.layers
    if args.bore > 0:                                  # destrava o bore: r1_min e bore_min seguem o requisito
        BORE = args.bore; args.bore_min = args.bore
    else:
        BORE = 0.090                                   # legado: r1_min=0.045+DIM/2+0.001 (bore ~92mm)
    ZPOS = np.arange(args.nz) * 0.020; NZ = args.nz   # comprimento axial configuravel
    print(f"[GA-homog] cubo={DIM*1e3:.0f}mm L={N_LAYERS}cam  objetivo=residuo shimmavel(SH<=3)  "
          f"B0 in [{args.b0_lo},{args.b0_hi}]T  ROI r<={ROI_R*1e3:.0f}mm  pop={args.pop} gen={args.gen} restarts={args.restarts}")
    champ = None; champ_score = np.inf
    for s in range(args.restarts):
        args.seed = s
        best, ev, cmask, Pproj = run(args)
        Bc, raw, uns, mass, bore, ni, ov = evaluate(best, ev, cmask, Pproj)
        feasible = (Bc >= args.b0_lo) and (bore >= args.bore_min) and (ov == 0)
        # score: residuo inalcancavel; descarta inviaveis com penalidade
        # selecao entre restarts: campo-aware (consistente com a fitness: premia B0 alto) + viabilidade
        score = uns/1000.0 + W_FIELD*(args.b0_hi - min(Bc, args.b0_hi))/args.b0_hi + (0 if feasible else 1e6)
        print(f"  [seed {s}] unshim={uns:.0f} raw={raw:.0f} B0={Bc*1e3:.1f}mT bore={bore*1e3:.1f}mm "
              f"N={ni} {'OK' if feasible else 'INVIAVEL'}")
        if score < champ_score:
            champ_score, champ = score, (best, Bc, raw, uns, mass, bore, ni)
    best, Bc, raw, uns, mass, bore, ni = champ
    radii, counts, grades = decode(best)
    print(f"[GA-homog] ESCOLHIDO: unshim={uns:.0f} raw={raw:.0f}ppm  B0={Bc*1e3:.1f}mT (f0={42.577*Bc:.3f}MHz)  "
          f"massa={mass:.1f}kg bore={bore*1e3:.1f}mm N={ni}  raios={np.round(radii*1e3,1)}mm")
    json.dump(to_json(best), open(args.o, "w"), indent=2, ensure_ascii=False)
    print(f"[GA-homog] escrito: {args.o}")
