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
// FEA colapso (fea_collapse_vessels.py + CalculiX): P_projeto=1,5x689=1033,5 bar (norma); escoamento governa.
// OD326 (pior caso) exigia t_estrut 18,1mm +CA3 = 21,1mm -> DECISAO: parede 22mm (=parede do magneto p/ ambos).
can_wall = 22; can_OD = mag_OD + 8 + 2*can_wall; can_L = mag_L + 40;   // can_OD=328 (cavidade mag_OD+8=284)

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

// ---- canister da ELETRONICA — em PÉ, base no PLATE, AO LADO do magneto (1-atm, <=40°C, RFPA+RFSoC+GPA+fonte) ----
// PROJETO_CANISTER_ELETRONICA §5.2: 1-atm Ti grau 5; P_projeto 1,5x689=1033,5 bar; CalculiX valida (sigma 539MPa<sy).
// FEA: OD200 exige so 14,1mm; mantido 16mm (mass-opt, SF_serv 2,22 >=1,5 por norma). (magneto e que foi p/ 22mm)
elec_OD = 200; elec_L = 350; elec_wall = 16;   // 16mm (mass-otimizado; CalculiX sigma 539MPa < sy 880)
elec_y  = can_OD/2 + elec_OD/2 + 25;     // ao lado do magneto (+Y), sobre o plate

// ---- DOCKING BASE + GUIA + STAB PLATE — conexao SOMENTE pela base (padrao SCM/SCMMB, API 17F) ----
// Modulo recuperavel por ROV pousa numa base; guide post central (alinhamento GROSSO=posicao) +
// 2 guide pins (alinhamento FINO/anti-rotacao); a stab plate (wet-mate power/sinal/fibra, familia
// Teledyne ODI Nautilus ROV/stab-mate) faz o MATE no ultimo curso; guide post = mandril de trava.
// As conexoes magneto<->eletronica saem pela BASE e correm POR BAIXO do skid ate o penetrador.
dock_OD = elec_OD + 40;  dock_H = 60;     // 240 — receiver baseplate sob o canister
post_D  = 50;  post_H = 150;              // guide post central (ponta conica) — alinh. grosso + mandril/trava
pin_D   = 25;  pin_H = 95;  pin_pcd = 170;// 2 guide pins anti-rotacao (alinh. fino), entram em buchas no modulo
stab_OD = 150; stab_t = 22;               // stab plate (placa de conexao wet-mate): macho no modulo, femea na base
explode = 120;                            // separacao p/ VISUALIZAR o acoplamento (=0 -> assentado p/ STEP)

dock_top = plate_top + dock_H;            // face superior da docking base
elec_zbase = dock_top + 120;              // canister pousa ELEVADO sobre a guia/stab (nao flange-flat)

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
module plate() color("slategray") {
    y0 = -mag_OD/2 - 80;                 // lado do magneto
    y1 = elec_y + dock_OD/2 + 10;        // estende ate cobrir a docking base (+Y)
    translate([0, (y0+y1)/2, plate_top - plate_t/2])
        cube([spool_L, y1-y0, plate_t], center=true);
}

// ----- berços (saddles) que sobem do plate e cradlam o magneto -----
// bases de suporte do magneto: 2 pedestais CILINDRICOS verticais; topo cortado pela carcaça (cradle = cilindro invertido); base rente ao plate (nada abaixo)
sup_D = 200; sup_x = can_L*0.42;
module saddles() color("darkkhaki")
    for (sx=[-1,1]) translate([sx*sup_x, 0, 0])
        difference() {
            translate([0,0, plate_top]) cylinder(h=-plate_top, d=sup_D);    // pedestal do plate (z=plate_top) ao eixo do magneto (z=0)
            rotate([0,90,0]) cylinder(h=sup_D+10, d=can_OD+1, center=true); // cradle = carcaça do magneto (cilindro invertido corta o topo)
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
module elec_canister()
    translate([0, elec_y, elec_zbase + elec_L/2 + explode]) {
        color("seagreen") difference() { cylinder(h=elec_L, d=elec_OD, center=true); cylinder(h=elec_L-2*elec_wall, d=elec_OD-2*elec_wall, center=true);}  // vaso 1-atm (colapso 689bar)
        // cold-plate de Cu/Al acoplado a parede (caminho termico ao mar) — fontes RFPA/GPA/DC-DC montam aqui
        color([0.72,0.45,0.20]) translate([0,0,0]) difference(){ cylinder(h=elec_L-2*elec_wall-10, d=elec_OD-2*elec_wall-2, center=true); cylinder(h=elec_L, d=elec_OD-2*elec_wall-14, center=true);}
        color("teal") translate([0,0,-elec_L/2-15]) cylinder(h=15, d=elec_OD*0.8, center=false);   // flange de base (interface da docking)
        color("mediumaquamarine") translate([0,0,elec_L/2]) cylinder(h=45, d1=elec_OD*0.55, d2=elec_OD*0.85, center=false);     // funil de pouso ROV (topo)
        color("gold") translate([0,0,elec_L/2+45]) rotate([0,90,0]) cylinder(h=elec_OD*0.6, d=18, center=true);      // handle ROV
        // ===== INTERFACE DE CONEXAO SO PELA BASE (desce/acopla na docking base) =====
        zb = -elec_L/2 - 15;                                                       // face inferior do flange
        color("mediumspringgreen") translate([0,0, zb-90]) difference(){           // funil de re-entry (captura o post)
            cylinder(h=90, d1=post_D+72, d2=post_D+14);
            translate([0,0,-1]) cylinder(h=92, d1=post_D+56, d2=post_D+4); }
        color("deepskyblue") translate([0,0, zb-6-stab_t/2]) difference(){         // stab plate MACHO (no modulo)
            cylinder(h=stab_t, d=stab_OD, center=true);
            cylinder(h=stab_t+2, d=post_D+10, center=true); }
        for (i=[0:3]) rotate([0,0,90*i]) color("dodgerblue")                        // pinos wet-mate (power/sinal/fibra)
            translate([stab_OD*0.3,0, zb-6-stab_t-14]) cylinder(h=18, d=12);
        for (sx=[-1,1]) color("khaki") translate([sx*pin_pcd/2, 0, zb-55]) difference(){ // buchas dos guide pins
            cylinder(h=58, d=pin_D+16); translate([0,0,-1]) cylinder(h=60, d=pin_D+4); }
    }

// ----- DOCKING BASE (receiver baseplate fixa no skid; equivalente ao SCMMB) -----
module docking_base() color("dimgray")
    translate([0, elec_y, plate_top + dock_H/2]) cylinder(h=dock_H, d=dock_OD, center=true);

// ----- GUIDE POST central (alinhamento grosso) + ponta conica + mandril de trava -----
module guide_post() color("orangered")
    translate([0, elec_y, dock_top]) {
        cylinder(h=post_H*0.72, d=post_D);
        translate([0,0,post_H*0.72]) cylinder(h=post_H*0.28, d1=post_D, d2=post_D*0.4);   // ponta conica
    }

// ----- 2 GUIDE PINS (alinhamento fino / anti-rotacao) -----
module guide_pins() color("yellow")
    for (sx=[-1,1]) translate([sx*pin_pcd/2, elec_y, dock_top]) cylinder(h=pin_H, d1=pin_D, d2=pin_D*0.55);

// ----- STAB PLATE FEMEA (receptaculo na base): furo central p/ post + 4 receptaculos wet-mate -----
module stab_female() color("blueviolet")
    translate([0, elec_y, dock_top + stab_t/2]) difference(){
        cylinder(h=stab_t, d=stab_OD, center=true);
        cylinder(h=stab_t+2, d=post_D+10, center=true);                                                  // passagem do post
        for (i=[0:3]) rotate([0,0,90*i]) translate([stab_OD*0.3,0,0]) cylinder(h=stab_t+2, d=16, center=true); // receptaculos
    }

// ----- conexoes pela BASE: BLOCO/raceway que CONTORNA (saddle) a carcaça do magneto -> docking base + penetrador + flying lead -----
chase_top = dock_top;        // topo do bloco RENTE a parte cinza (topo da docking base) — nao invade o flange roxo/stab
module connections() {
    yb0 = 0; yb1 = elec_y;   // da carcaça do magneto ate a docking base
    // BLOCO/raceway com SADDLE: contorna o shape cilindrico da carcaça (subtrai can_OD) ate a docking base
    color("olivedrab") difference() {
        translate([0, (yb0+yb1)/2, (plate_top + chase_top)/2])
            cube([100, yb1-yb0, chase_top-plate_top], center=true);
        rotate([0,90,0]) cylinder(h=130, d=can_OD+0.5, center=true);   // SADDLE = OD da carcaça do magneto (contorno colado)
    }
    // feedthrough da harness p/ o magneto — DENTRO do contorno do saddle (NAO exposto ao externo)
    color("royalblue") translate([0, 30, -165]) cylinder(h=30, d=34, center=true);
    color("navy")      translate([0, elec_y, plate_top - 20]) rotate([0,90,0]) cylinder(h=fr_L*0.8, d=22, center=true); // flying lead -> SCM/topside (sai da docking base)
}

// ===== assembly =====
flow_module();    // magneto+spool+pilha(RF/escudo/gradiente) DEITADO (horizontal)
plate();          // PLATE de base abaixo do magneto
saddles();        // berços do magneto sobre o plate
frame();          // skid no leito marinho (mudmat+postes+vigas+anodos+icamento)
elec_canister();  // eletronica EM PÉ, AO LADO do magneto, conexao SO pela base (recuperavel ROV)
docking_base();   // receiver baseplate sob o canister (SCMMB)
guide_post();     // guide post central (alinhamento grosso + mandril/trava)
guide_pins();     // 2 guide pins (alinhamento fino / anti-rotacao)
stab_female();    // stab plate femea (receptaculo wet-mate na base)
connections();    // harness por BAIXO do skid + penetrador de base + flying lead
