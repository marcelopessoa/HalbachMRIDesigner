#!/usr/bin/env python3
"""Mapas radiais de homogeneidade (plano transversal z=0) dos designs Halbach.
Desvio de |B0| em ppm na ROI do bore, via HalbachCylinder.calculateB (codigo real)."""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from HalbachCylinder import HalbachCylinder

DESIGNS = [
    ("examples/mpfm_03T_homog.json",      "336mT 9-slice (antigo)"),
    ("examples/mpfm_036T_2cam_13sl.json", "386mT 13-slice"),
    ("examples/mpfm_07T_champion.json",   "703mT 13-slice (ESCOLHIDO)"),
    ("examples/mpfm_036T_500mm.json",     "386mT 23-slice (500mm)"),
    ("examples/mpfm_07T_500mm.json",      "703mT 23-slice (500mm)"),
]
R = 0.020; N = 90

def fieldmap(js):
    cyl = HalbachCylinder(); cyl.loadJSON(js)
    xs = np.linspace(-R, R, N)
    XX, YY = np.meshgrid(xs, xs, indexing="ij")
    g = [XX.ravel()+1e-6, YY.ravel(), np.zeros(XX.size)]
    B = cyl.calculateB(g)
    mag = np.linalg.norm(B[:, :2], axis=1).reshape(N, N)
    mask = XX**2+YY**2 <= R*R
    bm = mag[mask].mean()
    ppm = (mag-bm)/bm*1e6; ppm[~mask] = np.nan
    return XX*1e3, YY*1e3, ppm, bm, np.nanmax(ppm)-np.nanmin(ppm)

def main():
    n = len(DESIGNS)
    fig, axs = plt.subplots(1, n, figsize=(4*n, 4.2))
    for ax, (js, ttl) in zip(axs, DESIGNS):
        try:
            X, Y, ppm, bm, pp = fieldmap(js)
        except Exception as e:
            ax.set_title(f"{ttl}\n(erro)"); continue
        v = np.nanmax(np.abs(ppm))
        c = ax.contourf(X, Y, ppm, levels=30, cmap="RdBu_r", vmin=-v, vmax=v)
        ax.contour(X, Y, ppm, levels=[-500,-100,0,100,500], colors="k", linewidths=0.4, alpha=0.5)
        th = np.linspace(0, 2*np.pi, 100)
        for rr in (10, 15, 20):
            ax.plot(rr*np.cos(th), rr*np.sin(th), "k--", lw=0.6, alpha=0.6)
        ax.set_aspect("equal"); ax.set_title(f"{ttl}\nB0={bm*1e3:.0f}mT  p-p={pp:.0f}ppm", fontsize=9)
        ax.set_xlabel("x [mm]");
        plt.colorbar(c, ax=ax, fraction=0.046, label="ppm")
    fig.suptitle("Mapas radiais de homogeneidade B0 (plano z=0, bruto pre-shim) — circulos r=10/15/20mm", fontsize=11)
    fig.tight_layout()
    fig.savefig("renders/homog_maps.png", dpi=110)
    print("renders/homog_maps.png")

if __name__ == "__main__":
    main()
