#!/usr/bin/env python3
"""Mapas radiais de homogeneidade POS-SHIM (plano z=0): aplica a cascata passivo+ativo
(shim_active.py) ao campo de cada design e plota o residuo em ppm."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from HalbachCylinder import HalbachCylinder
from shim_active import passive_stage, sensitivity_matrix, solve_currents

DESIGNS = [
    ("examples/mpfm_03T_homog.json",      "336mT 9-slice (antigo)"),
    ("examples/mpfm_036T_2cam_13sl.json", "386mT 13-slice"),
    ("examples/mpfm_07T_champion.json",   "703mT 13-slice (ESCOLHIDO)"),
    ("examples/mpfm_036T_500mm.json",     "386mT 23-slice (500mm)"),
    ("examples/mpfm_07T_500mm.json",      "703mT 23-slice (500mm)"),
]
R = 0.020; N = 70

def shimmed_map(js):
    cyl = HalbachCylinder(); cyl.loadJSON(js)
    xs = np.linspace(-R, R, N)
    XX, YY = np.meshgrid(xs, xs, indexing="ij")
    mask = (XX**2+YY**2 <= R*R)
    pts = np.stack([XX[mask]+1e-6, YY[mask], np.zeros(mask.sum())], axis=1)
    B = cyl.calculateB([pts[:, 0].copy(), pts[:, 1].copy(), pts[:, 2].copy()])
    b0v = B[:, :2].mean(axis=0); b0h = np.array([b0v[0], b0v[1], 0.0]); b0h /= np.linalg.norm(b0h)
    mag = B[:, :2] @ b0h[:2]; B0 = mag.mean()
    dppm = (mag-B0)/B0*1e6
    res_pass = passive_stage(pts, dppm, b0h)               # passivo (ferro)
    S_T = sensitivity_matrix(pts, b0h)
    _, res, _ = solve_currents(S_T, res_pass, B0)          # + ativo (multi-coil)
    grid = np.full((N, N), np.nan); grid[mask] = res
    return XX*1e3, YY*1e3, grid, B0, np.nanmax(grid)-np.nanmin(grid), np.sqrt(np.nanmean(res**2))

def main():
    n = len(DESIGNS)
    fig, axs = plt.subplots(1, n, figsize=(4*n, 4.3))
    for ax, (js, ttl) in zip(axs, DESIGNS):
        try:
            X, Y, g, B0, pp, rms = shimmed_map(js)
        except Exception as e:
            ax.set_title(f"{ttl}\n(erro {e})"); continue
        v = np.nanmax(np.abs(g)) or 1
        c = ax.contourf(X, Y, g, levels=30, cmap="RdBu_r", vmin=-v, vmax=v)
        th = np.linspace(0, 2*np.pi, 100)
        for rr in (10, 15, 20):
            ax.plot(rr*np.cos(th), rr*np.sin(th), "k--", lw=0.6, alpha=0.6)
        ax.set_aspect("equal"); ax.set_xlabel("x [mm]")
        ax.set_title(f"{ttl}\nB0={B0*1e3:.0f}mT  RMS={rms:.0f}  p-p={pp:.0f}ppm", fontsize=9)
        plt.colorbar(c, ax=ax, fraction=0.046, label="ppm")
    fig.suptitle("Homogeneidade B0 POS-SHIM (passivo ferro + ativo 32ch, plano z=0) — circulos r=10/15/20mm", fontsize=11)
    fig.tight_layout(); fig.savefig("renders/homog_maps_postshim.png", dpi=110)
    print("renders/homog_maps_postshim.png")

if __name__ == "__main__":
    main()
