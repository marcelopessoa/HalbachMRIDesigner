// ============================================================================
// coils_rf_gradient.scad — CAD conceitual das BOBINAS RF + GRADIENTE (DN65 bore104)
// RF: solenoide (B1 axial). Gradiente Gz: par de Maxwell (encode velocidade).
// Gradiente Gx: sela (Golay) p/ readout transversal. Geometria do anel concentrico:
//   RF a≈43,5mm (D87) · gradiente a_g≈47,5mm (D95) · escudo Faraday housing b≈50,5mm (D101).
// Refs: hw-components/RF (OCRA RF-coil), hw-components/Gradiente/GCoils (OCRA Rev005 X/Y/Z).
// ============================================================================
$fn = 64;
wire = 1.2;                 // raio do fio (representativo)

// --- RF solenoide ---
a_rf = 43.5; L_rf = 240; N_rf = 10;
// --- Gradiente Gz: par de Maxwell ---
a_g = 47.5; M_g = 6; sep = sqrt(3)*a_g;     // separacao otima de linearidade
// --- Gradiente Gx: sela (2 janelas, 120 graus) ---
a_sd = 47.5; z_sd = 70; arc = 120;

module torus(R, r) rotate_extrude() translate([R,0]) circle(r);

module solenoid(a, L, N) color([0.85,0.55,0.25])    // cobre
    for (i=[0:N-1]) translate([0,0,-L/2 + i*L/(N-1)]) torus(a, wire);

module maxwell(a, M, s) color("mediumpurple")        // Gz: 2 bundles, correntes opostas
    for (zc=[-s/2, s/2]) for (i=[0:M-1])
        translate([0,0, zc + (i-(M-1)/2)*2.6*wire]) torus(a, wire);

module saddle_window(a, zc, span)                    // 1 janela de sela: 2 arcos + 2 retas axiais
{
    // arcos (rotate_extrude com angulo) em z = +-zc
    for (zz=[-zc, zc]) translate([0,0,zz]) rotate_extrude(angle=span) translate([a,0]) circle(wire);
    // pernas axiais ligando os arcos, nas bordas do span
    for (ang=[0, span]) rotate([0,0,ang]) translate([a,0,0]) cylinder(h=2*zc, r=wire, center=true);
}
module saddle_Gx(a, zc, span) color("seagreen")      // Gx: 2 janelas opostas (+x e -x)
    for (rot=[0, 180]) rotate([0,0,rot - span/2]) saddle_window(a, zc, span);

// referencias de raio (translucidos)
module ref_tube(d_in, d_out, L, col)
    color(col, 0.12) difference(){ cylinder(h=L, d=d_out, center=true); cylinder(h=L+2, d=d_in, center=true); }

// ===== assembly das bobinas =====
solenoid(a_rf, L_rf, N_rf);             // RF (laranja/cobre)
maxwell(a_g, M_g, sep);                 // Gz (roxo)
saddle_Gx(a_sd, z_sd, arc);             // Gx (verde)
// contornos de referencia: liner Al2O3, escudo Faraday, bore
ref_tube(65, 84, L_rf, "ivory");        // liner Al2O3 (D65-84)
ref_tube(100, 102, L_rf, "orange");     // escudo Faraday housing (D100-102)
ref_tube(104, 105, L_rf, "gray");       // bore do magneto (D104)
