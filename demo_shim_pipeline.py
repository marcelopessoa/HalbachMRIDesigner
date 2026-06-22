#!/usr/bin/env python3
"""Demonstra o passo de shimming (§2/§5): extrai mapa B0 do design via codigo REAL
do repo (calculateB), monta mapa em ppm no plano transversal z=0 (ROI do tubo) e
resolve correcao por harmonicos esfericos com validation/shim_optimizer.py."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "validation"))
from HalbachCylinder import HalbachCylinder
from shim_optimizer import get_shim_basis, solve_shimming, calculate_ppm_metrics

def main(jsonfile, order=3, r_fov=0.020, n=64):
    cyl = HalbachCylinder()
    cyl.loadJSON(jsonfile)
    # plano transversal z=0
    xs = np.linspace(-r_fov, r_fov, n)
    XX, YY = np.meshgrid(xs, xs, indexing="ij")
    grid = [XX.ravel() + 1e-6, YY.ravel(), np.zeros(XX.size)]
    B = cyl.calculateB(grid)
    mag = np.linalg.norm(B[:, :2], axis=1).reshape(n, n)
    mask = (XX**2 + YY**2) <= r_fov**2
    Bmean = mag[mask].mean()
    fieldmap_ppm = (mag - Bmean) / Bmean * 1e6     # mapa em ppm
    fieldmap_ppm[~mask] = 0.0

    raw = calculate_ppm_metrics(fieldmap_ppm[mask])
    basis, names = get_shim_basis(shape=(n, n), order=order)
    coeffs, shimmed = solve_shimming(fieldmap_ppm, basis, mask=mask)
    sh = calculate_ppm_metrics(shimmed[mask])

    print(f"[shim] {os.path.basename(jsonfile)}  B0={Bmean*1000:.1f}mT  ROI disco r<=20mm  modos={names}")
    print(f"  BRUTO   : p-p={raw['peak_to_peak_ppm']:8.0f} ppm   rms={raw['rms_ppm']:7.0f} ppm")
    print(f"  SHIMMED : p-p={sh['peak_to_peak_ppm']:8.0f} ppm   rms={sh['rms_ppm']:7.0f} ppm   "
          f"(ganho {raw['peak_to_peak_ppm']/sh['peak_to_peak_ppm']:.1f}x)")
    print(f"  coef SH : " + ", ".join(f"{nm}={c:+.0f}" for nm, c in zip(names, coeffs)))

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "examples/mpfm_03T_opt.json",
         order=int(sys.argv[2]) if len(sys.argv) > 2 else 3)
