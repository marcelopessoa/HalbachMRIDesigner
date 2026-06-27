STEP gerados de halbach_assembly.scad (rev. colapso 1,5x: canister magneto parede 22mm/OD328, eletronica 16mm).
Fluxo: OpenSCAD .scad -> .csg ; FreeCAD 1.0 importCSG -> Import.export .step (BREP, 122 solidos).
- halbach_assembly.step : montagem completa (spool DN65, liner Al2O3, RF/gradiente/Faraday, canister magneto 1-atm, canister eletronica, frame/plate).
- halbach_assembly.csg  : intermediario.
coils_rf_gradient.scad NAO exportado (rotate_extrude do solenoide/sela falha no importCSG do FreeCAD).

coils_rf_gradient.step (tesselado, fallback STL->STEP): 51128 faces, fn=28.
 Fluxo: OpenSCAD .scad -> .stl ; FreeCAD Mesh.makeShapeFromMesh(tol 0.10mm) -> exportStep.
 BREP de malha (1 face por triangulo) -> arquivo grande (~55MB); use para visualizacao/montagem, nao para edicao parametrica.
