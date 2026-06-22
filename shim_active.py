#!/usr/bin/env python3
"""
shim_active.py — shim ATIVO multi-coil para o Halbach 0,30 T, acoplado ao rev_d_shim.

rev_d_shim (LincolnCB, verificado pos-pull) = driver de array multi-coil: 8 boards,
cada uma um AD5676 (DAC 8 ch) -> 64 canais; calibracao usa 32 ch. Escala DAC
(convert_waveform.py): +-5A -> 16-bit SIGNED +-32767 (dac=round(amps*32767/5)).
Formato [n_slices, n_channels]. Aqui: 32 bobinas (4 axiais x 8 azimutais).

Fisica-chave: no Halbach o B0 e' TRANSVERSAL (no plano xy). Só a componente do
campo da bobina ao longo de B0 (~x_hat) desloca a ressonancia. Logo a sensibilidade
de cada bobina = projecao do seu campo (Biot-Savart) em b0_hat, em ppm/A.

Pipeline: mapa B0 (calculateB) -> matriz de sensibilidade S (Biot-Savart) ->
LSQ com limites de corrente -> correntes por canal -> palavras DAC AD5676 +
CSV no formato rev_d_shim. Reporta residuo e compara com teto de 2a ordem.
"""
import sys
import numpy as np
from scipy.optimize import lsq_linear
from HalbachCylinder import HalbachCylinder

MU0 = 4e-7 * np.pi

# ---- array multi-coil (rev_d_shim: 32 canais) ----
N_AX = 4            # aneis axiais de bobinas
N_AZ = 8            # bobinas por anel (azimutal)
N_CH = N_AX * N_AZ  # 32 canais
FORMER_R = 0.046    # raio do former de shim: SOBRE o tubo do bore, dentro dos imas (48mm)
COIL_R = 0.016      # raio de cada loop de shim (cabe no former)
N_TURNS = 50        # espiras por bobina (shim coils reais sao multi-espira)
AX_SPAN = 0.120     # extensao axial coberta pelos aneis (+-60mm)
I_MAX = 5.0         # full-scale por canal (A): convert_waveform.py usa +-5A -> +-32767
N_SEG = 48          # segmentos por loop p/ Biot-Savart
DAC_FS = 32767      # AD5676 16-bit SIGNED (rev_d_shim: dac = round(amps*32767/5))


def coil_centers():
    """Centros e normais (radiais) das 32 bobinas no former."""
    zs = np.linspace(-AX_SPAN/2, AX_SPAN/2, N_AX)
    centers, normals, tang_u, tang_v = [], [], [], []
    for z in zs:
        for k in range(N_AZ):
            phi = 2*np.pi*k/N_AZ
            c = np.array([FORMER_R*np.cos(phi), FORMER_R*np.sin(phi), z])
            n = np.array([np.cos(phi), np.sin(phi), 0.0])      # normal radial
            u = np.array([0.0, 0.0, 1.0])                       # tangente axial
            v = np.cross(n, u)                                  # tangente azimutal
            centers.append(c); normals.append(n); tang_u.append(u); tang_v.append(v)
    return map(np.array, (centers, normals, tang_u, tang_v))


def loop_field(center, u, v, pts, n_seg=N_SEG):
    """Biot-Savart de um loop circular (raio COIL_R, 1 A) no plano (u,v) centrado
    em 'center'. Retorna B (n_pts,3) em Tesla."""
    th = np.linspace(0, 2*np.pi, n_seg, endpoint=False)
    dth = 2*np.pi/n_seg
    # pontos e tangentes do loop
    lp = center[None, :] + COIL_R*(np.cos(th)[:, None]*u[None, :] + np.sin(th)[:, None]*v[None, :])
    dl = COIL_R*dth*(-np.sin(th)[:, None]*u[None, :] + np.cos(th)[:, None]*v[None, :])
    B = np.zeros((pts.shape[0], 3))
    for s in range(n_seg):
        r = pts - lp[s][None, :]
        rn = np.linalg.norm(r, axis=1)
        rn3 = np.maximum(rn**3, 1e-12)
        B += MU0/(4*np.pi) * np.cross(dl[s][None, :], r) / rn3[:, None]
    return N_TURNS * B            # bobina multi-espira


# ---- estagio PASSIVO: ferro nas bandejas AXIAIS do former (r=47mm) ----
# Canais sao continuos em z -> ferro denso no eixo (carrega o perfil axial Z2/Z4/Z6,
# que o ativo de 4 aneis nao pega bem). Quantizacao fina = chapas finas empilhaveis.
PASS_R = 0.047; PASS_NAZ = 16; PASS_NAX = 9; PASS_AXSPAN = 0.150; PASS_QLEVELS = 12

def dipole_proj(center, m_hat, pts, b0_hat):
    """Campo (proj. em b0_hat) de um dipolo unitario em 'center' orientado m_hat.
    Unidade arbitraria (a amplitude vira o coeficiente do LSQ = quantidade de ferro)."""
    r = pts - center[None, :]
    rn = np.linalg.norm(r, axis=1); rn5 = np.maximum(rn**5, 1e-18); rn3 = np.maximum(rn**3, 1e-15)
    mdotr = r @ m_hat
    B = (3*mdotr[:, None]*r)/rn5[:, None] - m_hat[None, :]/rn3[:, None]
    return B @ b0_hat

def passive_stage(pts, dppm, b0_hat):
    """Emula shim passivo REALISTA: ferro so ADICIONA momento (susceptibilidade>0)
    ao longo de B0 -> coeficientes >=0 (lsq_linear bounded). Coarse de proposito;
    o ativo (signed, fino) limpa o residuo alcancavel. Retorna residuo (ppm)."""
    zs = np.linspace(-PASS_AXSPAN/2, PASS_AXSPAN/2, PASS_NAX)
    cols = []
    for z in zs:
        for k in range(PASS_NAZ):
            phi = 2*np.pi*k/PASS_NAZ
            c = np.array([PASS_R*np.cos(phi), PASS_R*np.sin(phi), z])
            cols.append(dipole_proj(c, b0_hat, pts, b0_hat))
    A = np.array(cols).T
    A = A/np.linalg.norm(A, axis=0)
    # ferro positivo + um offset global (frequencia/lock absorve o termo constante)
    A = np.hstack([A, np.ones((A.shape[0], 1))])
    nC = A.shape[1]
    # Tikhonov no ferro (evita superajuste do RMS que estragaria o resto); offset livre
    lam = 0.3
    Areg = np.vstack([A, np.sqrt(lam)*np.eye(nC)]); Areg[-1, -1] = 0
    breg = np.concatenate([-dppm, np.zeros(nC)])
    lb = np.concatenate([np.zeros(nC-1), [-np.inf]]); ub = np.full(nC, np.inf)
    coef = lsq_linear(Areg, breg, bounds=(lb, ub), max_iter=400).x
    # QUANTIZA o ferro (chapas finas empilhaveis ~PASS_QLEVELS niveis)
    qstep = coef[:-1].max()/PASS_QLEVELS if coef[:-1].max() > 0 else 1.0
    coef[:-1] = np.round(coef[:-1]/qstep)*qstep
    resid = dppm + A @ coef
    if rms(resid) >= rms(dppm):                # guarda em RMS (metrica do solver)
        return dppm
    return resid


def build_system(cyl, r_fov=0.020, z_fov=0.025, nr=7, nz=7):
    """ROI cilindrica do tubo. Retorna pontos, desvio de B0 em ppm, b0_hat, B0."""
    xs = np.linspace(-r_fov, r_fov, nr); zs = np.linspace(-z_fov, z_fov, nz)
    pts = np.array([(x, y, z) for x in xs for y in xs for z in zs
                    if x*x + y*y <= r_fov*r_fov])
    grid = [pts[:, 0].copy()+1e-6, pts[:, 1].copy(), pts[:, 2].copy()]
    B = cyl.calculateB(grid)
    b0_vec = B[:, :2].mean(axis=0)                       # direcao media de B0 (xy)
    b0_hat = np.array([b0_vec[0], b0_vec[1], 0.0]); b0_hat /= np.linalg.norm(b0_hat)
    mag = B[:, :2] @ b0_hat[:2]                          # componente ao longo de B0
    B0 = mag.mean()
    dppm = (mag - B0)/B0 * 1e6
    return pts, dppm, b0_hat, B0


def sensitivity_matrix(pts, b0_hat):
    """S[n_pts, N_CH]: ppm por Ampere de cada bobina (projecao em b0_hat)."""
    centers, normals, U, V = coil_centers()
    cols = []
    for i in range(N_CH):
        B = loop_field(centers[i], U[i], V[i], pts)
        cols.append(B @ b0_hat)        # T/A projetado
    return np.array(cols).T            # (n_pts, N_CH) em T/A


def solve_currents(S_T, dppm, B0, lam=0.5):
    """LSQ REGULARIZADO (Tikhonov) com limites: min ||S_ppm i + dppm||^2 + lam||i||^2,
    |i|<=I_MAX. A regularizacao evita overshoot contra modos inalcancaveis (que
    reduziriam RMS mas aumentariam o pico-a-pico). GUARDA: so aplica se reduzir p-p."""
    S_ppm = S_T / B0 * 1e6             # ppm/A
    n = S_ppm.shape[1]
    A = np.vstack([S_ppm, np.sqrt(lam)*np.eye(n)])
    b = np.concatenate([-dppm, np.zeros(n)])
    res = lsq_linear(A, b, bounds=(-I_MAX, I_MAX), max_iter=400)
    cur = res.x
    resid = dppm + S_ppm @ cur
    if rms(resid) >= rms(dppm):                # guarda em RMS (metrica do solver)
        cur = np.zeros(n); resid = dppm
    return cur, resid, S_ppm


def currents_to_dac(currents, i_max=I_MAX):
    """SIGNED 16-bit, identico a rev_d_shim convert_waveform.current_to_dac_value:
    dac = round(amps * 32767/5), clip +-32767."""
    code = np.round(currents * DAC_FS / i_max).astype(int)
    return np.clip(code, -DAC_FS, DAC_FS)


def pp(v): return v.max()-v.min()
def rms(v): return np.sqrt(np.mean(v**2))

def cascade(cyl, r_fov):
    pts, dppm, b0_hat, B0 = build_system(cyl, r_fov=r_fov)
    S_T = sensitivity_matrix(pts, b0_hat)
    res_pass = passive_stage(pts, dppm, b0_hat)
    cur, res, S_ppm = solve_currents(S_T, res_pass, B0)
    gam = 42.577e6                                  # Hz/T
    fwhm = 2.355 * gam * B0 * rms(res) * 1e-6       # FWHM da linha (Hz) ~ a partir do RMS
    t2s = 1/(np.pi*fwhm)*1e6 if fwhm > 0 else 0     # us
    return dict(pts=pts, B0=B0, dppm=dppm, res_pass=res_pass, res=res, cur=cur,
                S_ppm=S_ppm, fwhm=fwhm, t2s=t2s)

def main(jsonfile):
    cyl = HalbachCylinder(); cyl.loadJSON(jsonfile)
    pts0, _, b0h, B0 = build_system(cyl)
    auth = np.abs(sensitivity_matrix(pts0, b0h)/B0*1e6).mean(axis=0)
    print(f"[shim_active] {jsonfile.split('/')[-1]}  B0={B0*1e3:.1f}mT  "
          f"array {N_AX}x{N_AZ}={N_CH}ch @ r={FORMER_R*1e3:.0f}mm, {N_TURNS}esp, |i|<={I_MAX}A")
    print(f"  autoridade media {auth.mean():.0f} ppm/A (max {auth.max():.0f})")
    print(f"  CASCATA raw->passivo(ferro {PASS_NAZ}x{PASS_NAX}, {PASS_QLEVELS}niv)->ativo({N_CH}ch). METRICA=RMS (p-p secundario):\n")
    print(f"   {'ROI raio':>9} | {'RMS: raw':>8} {'passivo':>8} {'+ativo':>8} | {'(p-p ativo)':>11} {'FWHM':>8} {'T2*':>7} {'|i|max':>6}")
    last = None
    for r_fov in (0.020, 0.015, 0.010):
        c = cascade(cyl, r_fov)
        print(f"   r<={r_fov*1e3:4.0f}mm | {rms(c['dppm']):8.0f} {rms(c['res_pass']):8.0f} {rms(c['res']):8.0f} | "
              f"{pp(c['res']):11.0f} {c['fwhm']:6.0f}Hz {c['t2s']:5.0f}us {np.abs(c['cur']).max():5.2f}A")
        if r_fov == 0.015: last = c
    print("\n  (RMS em ppm = metrica otimizada/linewidth; o ativo+passivo reduzem RMS monotonico;")
    print("   p-p e' pior-caso na borda; alta-ordem azimutal decai p/ o centro -> melhor em ROI menor)\n")

    currents = last['cur']; dac = currents_to_dac(currents)
    np.savetxt("shim_active_amps.csv", currents[None, :], fmt="%.6f", delimiter=" ")
    np.savetxt("shim_active_dac.csv", dac[None, :], fmt="%d", delimiter=" ")
    print(f"  -> shim_active_amps.csv (Amps) + shim_active_dac.csv (AD5676 16-bit) [{N_CH} canais]")
    return currents, last['res']


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "examples/mpfm_03T_opt.json")
