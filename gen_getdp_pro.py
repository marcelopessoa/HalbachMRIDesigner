#!/usr/bin/env python3
"""Gera o .pro do GetDP a partir do .pickle de ImAs do HalbachMRIDesigner --fem
(o branch generateONELABfile do designer e' False; isto o substitui sem editar o
codigo vendorizado). Define angle_i/BR_i/mur_i + NumMagnets/SurfaceRegionOffset/
DSV/outputFilename e inclui o template magnetostatico MagSta_phi."""
import pickle
import sys

def main(project, dsv=0.04):
    with open(project + ".pickle", "rb") as f:
        s = pickle.load(f)
    dd = s["dataDict"]
    off = s["SurfaceRegionOffset"]
    lines = ["DefineConstant["]
    lines.append(f"  NumMagnets = {len(dd)},")
    lines.append(f"  SurfaceRegionOffset = {off},")
    lines.append(f"  DSV = {dsv},")
    lines.append(f'  outputFilename = "{project}",')
    body = []
    for d in dd:
        i = d["id"]
        body.append(f"  angle_{i} = {d['angle']}, BR_{i} = {d['BR']}, mur_{i} = {d['mur']}")
    lines.append(",\n".join(body))
    lines.append("];")
    lines.append('Include "templates/common_template.pro"')
    with open(project + ".pro", "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"[pro] {project}.pro escrito: {len(dd)} imas, DSV={dsv*1e3:.0f}mm")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "mpfm_03T_opt",
         float(sys.argv[2]) if len(sys.argv) > 2 else 0.04)
