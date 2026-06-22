import FreeCAD as App, importCSG, Import
BASE = r"\\wsl.localhost\Debian\home\marcelo\CS-MPFM-LFMRI\HalbachMRIDesigner"
importCSG.open(BASE + r"\halbach_assembly.csg"); doc=App.ActiveDocument; doc.recompute()
objs=[o for o in doc.Objects if hasattr(o,"Shape") and o.Shape.isValid() and o.Shape.Volume>1e-6]
Import.export(objs, BASE + r"\halbach_assembly.step")
App.Console.PrintMessage("\nASSEMBLY STEP OK solidos=%d\n"%len(objs))
