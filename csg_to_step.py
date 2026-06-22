import FreeCAD as App, importCSG, Import
BASE = r"\\wsl.localhost\Debian\home\marcelo\CS-MPFM-LFMRI\HalbachMRIDesigner"
importCSG.open(BASE + r"\mpfm_03T_fab.csg")
doc = App.ActiveDocument
doc.recompute()
objs = [o for o in doc.Objects if hasattr(o,"Shape") and o.Shape.Volume>1e-9 and o.Shape.isValid()]
App.Console.PrintMessage("\nObjetos solidos validos: %d\n" % len(objs))
Import.export(objs or doc.Objects, BASE + r"\mpfm_03T_fab.step")
App.Console.PrintMessage("STEP exportado\n")
