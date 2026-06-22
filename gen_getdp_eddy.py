#!/usr/bin/env python3
"""gen_getdp_eddy.py — SCAFFOLD do estudo de EDDY-CURRENT / atenuacao de B1 / perda (ΔQ) do
ESCUDO RF FENDIDO a f0 (~29,9 MHz), formulacao magnetodinamica time-harmonic (a*, GetDP).

Gera malha propria (NAO usa os imas): bobina RF (fonte) + escudo Cu fendido (condutor) + ar,
escreve o .pro a partir de templates/eddy_shield_template.pro e (opcional) roda o getdp.

CAVEAT (a validar): skin-depth Cu @29,9MHz ~12 um << parede real do escudo (18-35 um) e << tamanho
de elemento pratico. Aqui o escudo e' modelado como tubo de parede `--shield-wall` (default 0,5 mm,
sobre-espesso p/ meshabilidade) -> a PERDA absoluta nao e' quantitativa; serve p/ (a) validar o
pipeline e (b) comparar FENDIDO vs CONTINUO (a fenda e' o efeito de 1a ordem no B1). Para numero de
perda quantitativo: malha graduada na espessura real OU thin-shell/impedance BC (trabalho futuro).

Uso:  python gen_getdp_eddy.py [--slot/--no-slot] [--run]   (geometria do design DN65 bore104)
"""
import argparse, os, subprocess, sys
import numpy as np
import gmsh

# --- geometria (m), do halbach_assembly.scad (DN65 bore 104) ---
R_COIL_IN, R_COIL_OUT = 0.0425, 0.0445      # bobina RF (D85-89)
R_SH_IN                = 0.045               # escudo RF interno (D90)
COIL_L                 = 0.24                # comprimento dos enrolamentos
R_BOX                  = 0.12                # caixa de ar (semilado)
TAGS = dict(SHIELD=1, COIL=2, AIR=3, INF=4)


def build(slot=True, shield_wall=0.001, slot_w=0.003, msize=0.012):
    gmsh.initialize()
    gmsh.model.add("eddy")
    occ = gmsh.model.occ
    r_sh_out = R_SH_IN + shield_wall
    # tubos (cilindro externo - interno) via cut
    def tube(r_in, r_out, L):
        o = occ.addCylinder(0, 0, -L/2, 0, 0, L, r_out)
        i = occ.addCylinder(0, 0, -L/2, 0, 0, L, r_in)
        out, _ = occ.cut([(3, o)], [(3, i)])
        return out[0][1]
    coil = tube(R_COIL_IN, R_COIL_OUT, COIL_L)
    shield = tube(R_SH_IN, r_sh_out, COIL_L)
    if slot:                                  # fenda axial (paralela a B1): corta o tubo do escudo
        slotbox = occ.addBox(0, -slot_w/2, -COIL_L/2-0.001, r_sh_out*1.2, slot_w, COIL_L+0.002)
        out, _ = occ.cut([(3, shield)], [(3, slotbox)])
        shield = out[0][1]
    box = occ.addBox(-R_BOX, -R_BOX, -R_BOX, 2*R_BOX, 2*R_BOX, 2*R_BOX)
    occ.synchronize()
    # fragmenta tudo p/ interfaces conformes
    occ.fragment(occ.getEntities(3), [])
    occ.synchronize()
    # --- tagueamento ROBUSTO por RAIO DO BOUNDING-BOX (anel tem COM no eixo, mas bbox no raio externo) ---
    thr = 0.5*(R_COIL_OUT + R_SH_IN)          # limiar coil/shield (~midpoint do gap)
    shield_e, coil_e, air_e = [], [], []
    for (d, t) in occ.getEntities(3):
        bb = gmsh.model.getBoundingBox(3, t)
        rbb = max(abs(bb[0]), abs(bb[3]), abs(bb[1]), abs(bb[4]))   # max |x|,|y| do bbox
        if rbb > 0.5*R_BOX:                   # ocupa a caixa toda -> ar
            air_e.append(t)
        elif rbb > thr:                       # raio externo ~ escudo
            shield_e.append(t)
        else:                                 # raio externo ~ bobina
            coil_e.append(t)
    gmsh.model.addPhysicalGroup(3, shield_e, TAGS["SHIELD"])
    gmsh.model.addPhysicalGroup(3, coil_e, TAGS["COIL"])
    gmsh.model.addPhysicalGroup(3, air_e, TAGS["AIR"])
    # contorno externo (faces cujo centro esta na casca da caixa, em qq eixo = +-R_BOX)
    outer = []
    for (d, s) in gmsh.model.getBoundary([(3, t) for t in air_e], combined=True, oriented=False):
        cm = gmsh.model.occ.getCenterOfMass(2, abs(s))
        if any(abs(abs(cm[i]) - R_BOX) < 1e-3 for i in (0, 1, 2)):
            outer.append(abs(s))
    gmsh.model.addPhysicalGroup(2, outer, TAGS["INF"])
    print(f"[eddy-mesh] shield={len(shield_e)} coil={len(coil_e)} air={len(air_e)} outer={len(outer)} (slot={slot})")
    gmsh.option.setNumber("Mesh.MeshSizeMax", msize)
    # OOM-SAFE: NAO casar MeshSizeMin a shield_wall/2 (explode a malha -> 313k DOFs -> MUMPS OOM no WSL).
    # min sano (~2mm); o escudo fino fica sub-resolvido (ja era caveat) -> resultado QUALITATIVO.
    gmsh.option.setNumber("Mesh.MeshSizeMin", max(0.002, msize/4))
    gmsh.option.setNumber("Mesh.MshFileVersion", 2)
    gmsh.model.mesh.generate(3)
    gmsh.write("eddy_shield.msh")
    nn = gmsh.model.mesh.getNodes()[0].size
    gmsh.finalize()
    return nn


def write_pro(freq=29.9e6, sigma=5.96e7, js=1e6):
    lines = ["DefineConstant["]
    lines.append(f"  SHIELD={TAGS['SHIELD']}, COIL={TAGS['COIL']}, AIR={TAGS['AIR']}, INF={TAGS['INF']},")
    lines.append(f"  SIGMA={sigma}, FREQ={freq}, JS={js},")
    lines.append('  outputFilename = "eddy_shield"')
    lines.append("];")
    lines.append('Include "templates/eddy_shield_template.pro"')
    open("eddy_shield.pro", "w").write("\n".join(lines) + "\n")
    print("[eddy-pro] eddy_shield.pro escrito")


def _mem_limit(gb):
    """preexec p/ getdp: cap de memoria virtual (RLIMIT_AS) -> morre com ENOMEM em vez de derrubar o WSL."""
    import resource
    b = int(gb*1024**3)
    resource.setrlimit(resource.RLIMIT_AS, (b, b))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-slot", action="store_true", help="escudo CONTINUO (compara vs fendido)")
    ap.add_argument("--shield-wall", type=float, default=0.001)
    ap.add_argument("--msize", type=float, default=0.012)
    ap.add_argument("--max-nodes", type=int, default=40000, help="guarda OOM: recusa solve se a malha exceder")
    ap.add_argument("--mem-gb", type=float, default=6.0, help="cap de RAM do getdp (RLIMIT_AS)")
    ap.add_argument("--iter", action="store_true", help="solver ITERATIVO (gmres+ilu) em vez de MUMPS LU")
    ap.add_argument("--run", action="store_true", help="roda getdp apos gerar")
    a = ap.parse_args()
    nn = build(slot=not a.no_slot, shield_wall=a.shield_wall, msize=a.msize)
    write_pro()
    print(f"[eddy-mesh] nodes={nn} (guarda max-nodes={a.max_nodes})")
    if a.run:
        if nn > a.max_nodes:
            print(f"[eddy-run] ABORTADO: malha {nn} nós > --max-nodes {a.max_nodes}. "
                  f"Aumente --msize (mais grosseira) ou --max-nodes (e RAM). 3D edge eddy LU complexo "
                  f"escala mal -> p/ numero quantitativo use thin-shell/impedance BC (ver template/doc).")
            sys.exit(2)
        getdp = "../.venv-halbach/bin/getdp"
        getdp = getdp if os.path.exists(getdp) else "getdp"
        cmd = [getdp, "eddy_shield.pro", "-msh", "eddy_shield.msh", "-solve", "Eddy", "-pos", "map"]
        if a.iter:
            cmd += ["-ksp_type", "gmres", "-pc_type", "ilu", "-ksp_rtol", "1e-6", "-ksp_max_it", "2000"]
        print(f"[eddy-run] getdp (cap {a.mem_gb}GB, {'iterativo' if a.iter else 'MUMPS LU'}) ...")
        subprocess.run(cmd, check=False, preexec_fn=lambda: _mem_limit(a.mem_gb))
