#!/usr/bin/env python3
"""Valida a config MPFM usando o codigo REAL do repo (HalbachCylinder.loadJSON +
calculateB), mas avaliando na ROI correta do tubo (o HalbachMRIDesigner.py fixa
dsv=0.2 m, esfera grande demais p/ bore pequeno -> engloba os imas -> singularidade).
ROI: cilindro FOV radial +-20 mm, axial +-25 mm (regiao de escoamento util)."""
import sys
import numpy as np
from HalbachCylinder import HalbachCylinder

def main(jsonfile):
    cyl = HalbachCylinder()
    cyl.loadJSON(jsonfile)
    n_imas = sum(len(r.magnets) for s in cyl.slices for r in s.rings)
    print(f"[load] {jsonfile}: {len(cyl.slices)} slices, {n_imas} imas")

    # ROI do tubo: grade cilindrica
    r_fov, z_fov = 0.020, 0.025
    xs = np.linspace(-r_fov, r_fov, 9)
    zs = np.linspace(-z_fov, z_fov, 11)
    pts = [(x, y, z) for x in xs for y in xs for z in zs
           if x*x + y*y <= r_fov*r_fov]
    pts = np.array(pts, dtype=np.float64)
    # evita r=0 exato (formula dipolar) deslocando epsilon
    grid = [pts[:, 0].copy() + 1e-6, pts[:, 1].copy(), pts[:, 2].copy()]

    B = cyl.calculateB(grid)
    mag = np.linalg.norm(B[:, :2], axis=1)            # B0 transversal (Halbach)
    center = (np.abs(pts[:, 0]) < 1e-3) & (np.abs(pts[:, 1]) < 1e-3) & (np.abs(pts[:, 2]) < 1e-3)
    Bc = mag[center].mean() if center.any() else mag.mean()
    ppm = (mag.max() - mag.min()) / mag.mean() * 1e6
    f0 = 42.577e6 * Bc                                 # Larmor 1H

    print(f"[campo] ROI cilindro r<=20mm, z<=+-25mm ({pts.shape[0]} pts)")
    print(f"  B0 central      = {Bc*1000:8.1f} mT")
    print(f"  B0 medio        = {mag.mean()*1000:8.1f} mT")
    print(f"  f0(1H)          = {f0/1e6:8.3f} MHz")
    print(f"  homogeneidade   = {ppm:8.0f} ppm p-p (BRUTO, pre-shim)")
    print(f"  >> derating FEM (mur=1.05+demag, -10..20%) -> {Bc*0.85*1000:.0f}..{Bc*0.90*1000:.0f} mT")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "examples/mpfm_03T_dn50.json")
