#!/usr/bin/env python3
"""Shim ATIVO 3D — quantifica o residuo alcancavel com um jogo de 2a ordem
(8 termos = 8 canais DAC do rev_d_shim). Campo via codigo real do repo
(HalbachCylinder.calculateB) numa ROI cilindrica do tubo; ajuste por minimos
quadrados na base de harmonicos esfericos solidos reais ate ordem 2.
Z0 (offset) NAO entra: e' absorvido pela frequencia do sintetizador / lock."""
import sys
import numpy as np
from HalbachCylinder import HalbachCylinder

# Harmonicos esfericos solidos REAIS (cartesianos), ate ordem 2. Z0 excluido.
def sh_basis(x, y, z):
    return {
        "Z (lin)":   z,
        "X (lin)":   x,
        "Y (lin)":   y,
        "Z2":        2*z*z - (x*x + y*y),      # 2z^2 - x^2 - y^2
        "ZX":        x*z,
        "ZY":        y*z,
        "X2-Y2":     x*x - y*y,
        "XY":        x*y,
    }

def main(jsonfile, order_set="full"):
    cyl = HalbachCylinder(); cyl.loadJSON(jsonfile)
    # ROI cilindrica do tubo
    r_fov, z_fov = 0.020, 0.025
    xs = np.linspace(-r_fov, r_fov, 9); zs = np.linspace(-z_fov, z_fov, 9)
    P = np.array([(x, y, z) for x in xs for y in xs for z in zs
                  if x*x + y*y <= r_fov*r_fov])
    grid = [P[:, 0].copy() + 1e-6, P[:, 1].copy(), P[:, 2].copy()]
    B = cyl.calculateB(grid)
    mag = np.linalg.norm(B[:, :2], axis=1)
    B0 = mag.mean()
    dppm = (mag - B0) / B0 * 1e6            # desvio em ppm

    basis = sh_basis(P[:, 0], P[:, 1], P[:, 2])
    names = list(basis.keys())
    if order_set == "linear":
        names = names[:3]
    A = np.stack([basis[n] for n in names], axis=1)
    # normaliza colunas p/ condicionamento
    scale = np.linalg.norm(A, axis=0); A = A / scale
    coef, *_ = np.linalg.lstsq(A, -dppm, rcond=None)
    resid = dppm + A @ coef

    pp0 = dppm.max() - dppm.min(); rms0 = np.sqrt(np.mean(dppm**2))
    pp1 = resid.max() - resid.min(); rms1 = np.sqrt(np.mean(resid**2))
    print(f"[shim3D] {jsonfile.split('/')[-1]}  B0={B0*1e3:.1f}mT  ROI cil r<=20 z<=+-25mm ({P.shape[0]} pts)  termos={names}")
    print(f"  BRUTO    : p-p={pp0:8.0f} ppm   rms={rms0:7.0f} ppm")
    print(f"  SHIMMED  : p-p={pp1:8.0f} ppm   rms={rms1:7.0f} ppm   (ganho p-p {pp0/pp1:.1f}x, rms {rms0/rms1:.1f}x)")

if __name__ == "__main__":
    f = sys.argv[1] if len(sys.argv) > 1 else "examples/mpfm_03T_opt.json"
    main(f, sys.argv[2] if len(sys.argv) > 2 else "full")
