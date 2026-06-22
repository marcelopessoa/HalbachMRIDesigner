// ============================================================================
// Assembly conceitual do MPFM (Projeto Detalhado) — magneto 0,7T DEITADO (horizontal)
// sobre FRAME/skid no leito marinho; canister de eletronica PENDURADO no frame.
// Eixo de fluxo = X (horizontal). Z = vertical (cima). Geometria LIMPA -> STEP-able. mm.
// CONFORMIDADE: API 17S §5.7 (canister eletronica <=40°C, recuperavel), §5.11, API 6A (flanges).
// REPROJETO DN65 v2 (2026-06-22): flow bore 2-9/16" API 6A; bore magneto 104mm; barreira de pressao
// = LINER CERAMICO Al2O3 (NAO-condutor, RF-TRANSPARENTE; metal bloquearia o B1 -> skin-depth ~100um
// @29,9MHz). nz=9 (17 slices) p/ recuperar 0,7T. PILHA CONCENTRICA: RF + escudo RF (fenda ||B1) +
// gradiente + escudo gradiente (fenda ortogonal). Ver ../industrial-production/prototipo/REPROJETO_MAGNETO_BORE_2026.md
// ============================================================================

// ---- magneto (escolhido DN65 nz9: 3 cam x 17 slices, GA 0,704T / FEM 0,725T, bore 104mm) ----
mag_OD = 276; mag_L = 340; mag_bore = 104;     // nz=9 -> 17 slices -> mag_L 340; bore 104

// ---- LINER CERAMICO Al2O3 (689 bar, RF-TRANSPARENTE) + flow bore (horizontal) ----
flow_ID = 65; liner_wall = 9.5;         // DN65=2-9/16" API 6A; Al2O3 (FICHA: P_max>1000bar, RF-transparente)
liner_OD = flow_ID + 2*liner_wall;      // 84 (ceramica NAO-condutora: deixa o RF passar; metal bloquearia)
spool_L  = mag_L + 220;                 // estende p/ as flanges metalicas API 6A (axis X apos deitar)
roi_ID   = 30;                          // ROI homogenea r<=15mm = nucleo de medicao

// ---- PILHA CONCENTRICA no anel (anel = (104-84)/2 = 10mm/lado; pilha ~8,5mm + folga) ----
//   ordem do fluido p/ fora: liner Al2O3 | folga | RF coil | escudo RF | folga | gradiente | escudo grad | folga | bore
coil_L  = 240;                          // comprimento dos enrolamentos (dentro do bore)
slit_w  = 3;                            // largura da fenda anti-eddy (overlap capacitivo conceitual)
rf_in   = liner_OD + 1;  rf_out  = rf_in + 4;     // RF solenoide + former (D 85->89), a~43,5mm
gr_in   = rf_out + 3;    gr_out  = gr_in + 5;      // bobina de gradiente ENTRE coil e escudo (D 92->97)
fsh_in  = gr_out + 3;    fsh_out = fsh_in + 2;     // ESCUDO FARADAY = housing externo, 32 fendas axiais ||B1 (D 100->102, b~50,5; b/a~1,16; folga ate bore 104)

// ---- conexoes API 6A (flange 2-9/16" / 3-1/8" classe alta) ----
flange_OD = 210; flange_t = 44; bolt_pcd = 168; bolt_d = 24; n_bolts = 8;

// ---- carcaça/canister do magneto (1-atm) ----
can_OD = 326; can_wall = 14; can_L = mag_L + 40;

// ---- FRAME / skid no leito marinho ----
fr_L = spool_L + 120;     // comprimento (ao longo de X)
fr_W = 720;               // largura (Y) — acomoda magneto + canister eletronica ao lado
fr_H = 520;               // altura (Z) do mudmat ao topo
beam = 50;                // secao das vigas
mud_t = 30;               // espessura do mudmat

// ---- PLATE de base (abaixo do magneto; carrega magneto em berços + canister eletronica) ----
plate_t   = 30;
plate_top = -can_OD/2 - 40;          // logo abaixo do magneto
plate_y0  = -fr_W/2 + beam;          // estende ate o lado +Y onde fica a eletronica

// ---- canister da ELETRONICA — em PÉ, base apoiada no PLATE, AO LADO do magneto (1-atm, <=40°C) ----
elec_OD = 180; elec_L = 300;
elec_y  = can_OD/2 + elec_OD/2 + 25;     // ao lado do magneto (+Y), sobre o plate
elec_zbase = plate_top;                   // base apoiada no plate

$fn = 96;

// ----- modulo de fluxo (construido em torno do eixo Z; depois deitado p/ X) -----
module spool_z() {
    difference() {
        union() {
            cylinder(h=spool_L, d=liner_OD, center=true);
            for (z=[-1,1]) translate([0,0,z*(spool_L/2-flange_t/2)]) cylinder(h=flange_t, d=flange_OD, center=true);
        }
        cylinder(h=spool_L+2, d=flow_ID, center=true);
        for (z=[-1,1]) translate([0,0,z*(spool_L/2-flange_t/2)])
            for (i=[0:n_bolts-1]) rotate([0,0,i*360/n_bolts]) translate([bolt_pcd/2,0,0]) cylinder(h=flange_t+2, d=bolt_d, center=true);
    }
}
module canister_z() {
    difference() {
        cylinder(h=can_L, d=can_OD, center=true);
        cylinder(h=can_L-2*can_wall, d=mag_OD+8, center=true);
        cylinder(h=can_L+2, d=liner_OD+6, center=true);
    }
}
module magnet_z() { difference() { cylinder(h=mag_L, d=mag_OD, center=true); cylinder(h=mag_L+2, d=mag_bore, center=true);} }
module roi_z()    { difference() { cylinder(h=mag_L, d=roi_ID+4, center=true); cylinder(h=mag_L+2, d=roi_ID, center=true);} }

// ----- PILHA CONCENTRICA: RF + escudos fendidos + gradiente -----
module tube_z(d_in,d_out,L) difference(){ cylinder(h=L,d=d_out,center=true); cylinder(h=L+2,d=d_in,center=true); }
module rf_coil_z()   tube_z(rf_in, rf_out, coil_L);
module grad_coil_z() tube_z(gr_in, gr_out, coil_L);
// ESCUDO FARADAY = housing externo; 32 fendas axiais ||B1 (pontes capacitivas de overlap nao desenhadas)
module faraday_z()   difference(){ tube_z(fsh_in, fsh_out, coil_L);
    for (i=[0:31]) rotate([0,0,i*360/32]) translate([fsh_out/2,0,0]) cube([4, slit_w, coil_L+2], center=true); }

module flow_module() {           // DEITADO: eixo de fluxo ao longo de X
    rotate([0,90,0]) {
        color("steelblue")          canister_z();
        color([0.93,0.93,0.88])     spool_z();        // LINER Al2O3 ceramico (RF-transparente; flanges metalicas)
        color("gold")               magnet_z();
        color("tomato")             roi_z();
        // --- pilha concentrica no anel ---
        color([0.85,0.55,0.25])     rf_coil_z();       // RF (cobre), a~43,5mm
        color("mediumpurple")       grad_coil_z();      // gradiente (entre coil e escudo)
        color("orange")             faraday_z();        // ESCUDO FARADAY = housing externo, 32 fendas axiais ||B1 (b/a~1,16)
    }
}

// ----- PLATE de base abaixo do magneto (carrega magneto + canister eletronica) -----
module plate() color("slategray")
    translate([0, (plate_y0)/2, plate_top - plate_t/2])
        cube([spool_L, fr_W - 2*beam, plate_t], center=true);

// ----- berços (saddles) que sobem do plate e cradlam o magneto -----
module saddles() color("slategray")
    for (sx=[-1,1]) translate([sx*can_L*0.35, 0, (plate_top + (-can_OD/2))/2])
        difference() {
            cube([60, can_OD*0.9, abs(plate_top + can_OD/2)], center=true);
            translate([0,0,abs(plate_top + can_OD/2)/2]) rotate([0,90,0]) cylinder(h=70, d=can_OD+6, center=true);
        }

// ----- frame / skid no leito marinho -----
module frame() color("dimgray") {
    translate([0,0,-fr_H/2]) cube([fr_L, fr_W, mud_t], center=true);                                  // mudmat
    for (sx=[-1,1], sy=[-1,1]) translate([sx*(fr_L/2-beam/2), sy*(fr_W/2-beam/2), 0]) cube([beam,beam,fr_H], center=true); // postes
    for (sy=[-1,1], sz=[-1,1]) translate([0, sy*(fr_W/2-beam/2), sz*(fr_H/2-beam/2)]) cube([fr_L,beam,beam], center=true); // vigas long.
    for (sx=[-1,1]) translate([sx*(fr_L/2-beam/2),0,fr_H/2-beam/2]) cube([beam,fr_W,beam], center=true);                   // vigas transv. topo
    for (sx=[-1,1]) translate([sx*(fr_L/2-40),fr_W/2-30,-fr_H/2+mud_t/2+40]) cube([40,80,60], center=true);                // anodos
    for (sx=[-1,1]) translate([sx*(fr_L/2-beam),0,fr_H/2]) rotate([90,0,0]) cylinder(h=fr_W*0.5,d=24,center=true);         // icamento ROV
}

// ----- canister de eletronica EM PÉ, base no PLATE, ao lado do magneto; recuperavel por ROV -----
module elec_canister() color("seagreen")
    translate([0, elec_y, elec_zbase + elec_L/2]) {
        difference() { cylinder(h=elec_L, d=elec_OD, center=true); cylinder(h=elec_L-24, d=elec_OD-24, center=true);}
        translate([0,0,-elec_L/2-15]) cylinder(h=15, d=elec_OD*0.8, center=false);                    // flange de base no plate
        translate([0,0,elec_L/2]) cylinder(h=45, d1=elec_OD*0.55, d2=elec_OD*0.85, center=false);     // funil de pouso ROV (topo)
        translate([0,0,elec_L/2+45]) rotate([0,90,0]) cylinder(h=elec_OD*0.6, d=18, center=true);      // handle ROV
    }

// ----- conexoes: wet-mate eletrico+fibra (eletronica<->magneto) + penetrador + flying lead -----
module connections() color("orange") {
    translate([0, (elec_y + can_OD/2)/2, 0]) rotate([90,0,0]) cylinder(h=elec_y-can_OD/2, d=40, center=true); // jumper wet-mate
    translate([0, can_OD/2-5, 0]) rotate([90,0,0]) cylinder(h=30, d=36, center=true);                          // penetrador (sela 1-atm)
    translate([fr_L/2-30, 0, -fr_H/2+mud_t+30]) rotate([0,90,0]) cylinder(h=120, d=22, center=true);           // flying lead -> SCM/topside
}

// ===== assembly =====
flow_module();    // magneto+spool+pilha(RF/escudo/gradiente) DEITADO (horizontal)
plate();          // PLATE de base abaixo do magneto
saddles();        // berços do magneto sobre o plate
frame();          // skid no leito marinho (mudmat+postes+vigas+anodos+icamento)
elec_canister();  // eletronica EM PÉ, base no plate, AO LADO do magneto (recuperavel ROV)
connections();    // wet-mate + penetrador + flying lead
