// ============================================================================
// eddy_shield_template.pro — MAGNETODINAMICA time-harmonic (formulacao a*, 3D edge)
// Estuda eddy-current / atenuacao de B1 / perda (ΔQ) do ESCUDO RF fendido a f0 (~29,9 MHz).
// Gerado/parametrizado por gen_getdp_eddy.py (DefineConstant abaixo).
//
// SCAFFOLD (jun/2026) — A VALIDAR. Caveat fisico: skin-depth do Cu @29,9MHz ~12 um << parede do
// escudo (18-35 um) e << tamanho de elemento pratico -> resolver a pele em volume exige malha
// finissima OU condicao de contorno de impedancia (thin-shell). Este template resolve em VOLUME
// (massa condutora); para o escudo fino real, usar malha graduada na espessura ou thin-shell BC.
// Regioes (tags do gmsh, via gen_getdp_eddy.py):
//   SHIELD = escudo Cu (condutor, eddy)  |  COIL = bobina fonte (Js imposta)
//   AIR = ar  |  INF = contorno externo (a x n = 0, Dirichlet)
// ============================================================================
DefineConstant[
  SHIELD = 1, COIL = 2, AIR = 3, INF = 4,
  SIGMA = 5.96e7,     // condutividade do escudo (Cu)
  FREQ  = 29.9e6,     // frequencia de operacao (1H @ 0,7T)
  JS    = 1.0e6,      // densidade de corrente fonte na bobina (A/m^2, normalizada)
  outputFilename = "eddy_shield"
];

Group {
  Shield = Region[SHIELD];
  Coil   = Region[COIL];
  Air    = Region[AIR];
  Sur_Inf = Region[INF];
  Vol_C  = Region[{Shield}];                 // condutores (eddy currents)
  Vol_S  = Region[{Coil}];                   // fonte (corrente imposta)
  Vol_NonC = Region[{Air, Coil}];            // nao-condutores (precisam de gauge)
  Vol_Mag = Region[{Shield, Air, Coil}];     // dominio total
}

Function {
  mu0 = 4*Pi*1e-7;
  nu[] = 1.0/mu0;
  sigma[Shield] = SIGMA;
  // solenoide: corrente azimutal (no plano xy), magnitude JS
  js[Coil] = JS * Vector[ -Sin[Atan2[Y[],X[]]], Cos[Atan2[Y[],X[]]], 0. ];
}

Jacobian { { Name Vol; Case { { Region All; Jacobian Vol; } } } }
Integration { { Name I1; Case { { Type Gauss; Case {
  { GeoElement Tetrahedron; NumberOfPoints 4; }
  { GeoElement Triangle;    NumberOfPoints 3; } } } } } }

Constraint {
  { Name a_dirichlet; Case { { Region Sur_Inf; Value 0.; } } }
  // gauge de arvore (tree-cotree) nos nao-condutores: fixa as arestas da co-arvore
  { Name a_gauge; Case { { Region Vol_Mag; SubRegion Vol_S; } } }
}

FunctionSpace {
  { Name Hcurl_a; Type Form1;
    BasisFunction { { Name se; NameOfCoef ae; Function BF_Edge; Support Vol_Mag; Entity EdgesOf[All]; } }
    Constraint {
      { NameOfCoef ae; EntityType EdgesOf; NameOfConstraint a_dirichlet; }
      { NameOfCoef ae; EntityType EdgesOfTreeIn; EntitySubType StartingOn;
        NameOfConstraint a_gauge; }    // gauge: arvore comeca na fonte
    }
  }
}

Formulation {
  { Name MagDyn_a; Type FemEquation;
    Quantity { { Name a; Type Local; NameOfSpace Hcurl_a; } }
    Equation {
      Galerkin { [ nu[] * Dof{d a} , {d a} ]; In Vol_Mag; Jacobian Vol; Integration I1; }
      Galerkin { DtDof[ sigma[] * Dof{a} , {a} ]; In Vol_C; Jacobian Vol; Integration I1; }
      Galerkin { [ -js[] , {a} ]; In Vol_S; Jacobian Vol; Integration I1; }
    }
  }
}

Resolution {
  { Name Eddy;
    System { { Name A; NameOfFormulation MagDyn_a; Type Complex; Frequency FREQ; } }
    Operation { Generate[A]; Solve[A]; SaveSolution[A]; PostOperation[map]; }
  }
}

PostProcessing {
  { Name map; NameOfFormulation MagDyn_a;
    Quantity {
      { Name b;     Value { Term { [ {d a} ]; In Vol_Mag; Jacobian Vol; } } }
      { Name b_abs; Value { Term { [ Norm[{d a}] ]; In Vol_Mag; Jacobian Vol; } } }
      { Name jeddy; Value { Term { [ -sigma[]*Dt[{a}] ]; In Vol_C; Jacobian Vol; } } }
      // perda Joule total (-> ΔQ): P = 1/2 ∫ sigma |ω a|^2  (time-harmonic, media)
      { Name ploss; Value { Integral { [ 0.5*sigma[]*SquNorm[Dt[{a}]] ]; In Vol_C; Jacobian Vol; Integration I1; } } }
    }
  }
}

PostOperation {
  { Name map; NameOfPostProcessing map;
    Operation {
      Print[ b_abs, OnGrid {$B*Cos[$A], $B*Sin[$A], 0} {0:2*Pi:2*Pi/360, 0:0.015:0.015/30, 0},
             File StrCat[outputFilename, "_B1_frontal.pos"] ];
      Print[ ploss[Vol_C], OnGlobal, Format Table, File StrCat[outputFilename, "_ploss.txt"] ];
    }
  }
}
