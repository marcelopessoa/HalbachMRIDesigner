// ============================================================================
// Suporte / cradle imprimivel para o Halbach 0,30 T (CS-MPFM-LFMRI)
// Geometria LIMPA (so cylinder/cube/rotate/difference) -> converte p/ STEP no
// FreeCAD sem os linear_extrude problematicos do stand original.
// Inclui ESPACOS INTELIGENTES para SHIM PASSIVO (canais axiais p/ ferro/imas)
// conforme conceito Siemens ("shim iron em bandejas"). Imprimir no CIMATEC.
// Unidades: mm.
// ============================================================================

// ---- parametros do magneto (espelham examples/mpfm_03T_opt.json) ----
n_slices        = 9;
slice_sep       = 20;      // passo axial das slices
slice_thick     = 12;      // espessura de cada disco-holder
slice_OD        = 190;     // 2 x outerRadius (95 mm)
bore_free_r     = 45;      // raio livre do bore (tubo DN50 + RF + Faraday)

// ---- suporte ----
plate_t         = 10;      // espessura das placas de extremidade
rod_n           = 8;       // tirantes longitudinais
rod_r           = 6;       // raio dos tirantes
rod_pcd         = 100;     // raio do circulo de tirantes (PCD/2)
shell_clear     = 3;       // folga radial placa->slice
foot_h          = 30;      // altura do pe de bancada
foot_w          = 180;

// ---- SHIM PASSIVO: canais axiais (bandejas) ----
// Former fino entre o RF/Faraday (r<43) e a face interna dos imas (r~47,7 do design
// escolhido mpfm_03T_homog). 16 canais AZIMUTAIS continuos em z -> ferro denso no eixo
// (modelo de shim usa 16x9=144 posicoes: 9 niveis z por canal). Recuado p/ NAO colidir
// com os imas internos (antes r47+parede4=51 > 47,7 -> colisao; agora 43+3=46, livra).
shim_n          = 16;      // n. de bandejas azimutais (1 a cada 22,5 deg)
shim_former_r   = 43;      // raio interno do former de shim passivo
shim_slot_w     = 6;       // largura da bandeja
shim_slot_d     = 3;       // profundidade radial (aceita chapa de aco/ima fino)
shim_former_t   = 3;       // parede do former (externo em 46mm, livra imas em 47,7)

L = (n_slices-1)*slice_sep + slice_thick;   // comprimento util
plate_z = L/2 + plate_t/2;                  // posicao das placas

$fn = 120;

// ---------------------------------------------------------------------------
module end_plate() {
    difference() {
        cylinder(h=plate_t, r=slice_OD/2 + shell_clear + 8, center=true);
        // bore central
        cylinder(h=plate_t+2, r=bore_free_r, center=true);
        // alojamento dos tirantes
        for (i=[0:rod_n-1]) rotate([0,0,i*360/rod_n])
            translate([rod_pcd,0,0]) cylinder(h=plate_t+2, r=rod_r+0.3, center=true);
        // rebaixo de registro p/ a slice (encaixe)
        cylinder(h=plate_t, r=slice_OD/2 + 0.4, center=false);
    }
}

module tie_rods() {
    for (i=[0:rod_n-1]) rotate([0,0,i*360/rod_n])
        translate([rod_pcd,0,0]) cylinder(h=L+2*plate_t, r=rod_r, center=true);
}

// Former de shim passivo com bandejas axiais indexadas (espacos inteligentes)
module shim_former() {
    difference() {
        cylinder(h=L, r=shim_former_r+shim_former_t, center=true);
        cylinder(h=L+2, r=shim_former_r, center=true);             // vazio interno (bore)
        // bandejas axiais: canais retangulares passantes p/ inserir chapa/ima
        for (i=[0:shim_n-1]) rotate([0,0,i*360/shim_n])
            translate([shim_former_r+shim_former_t/2, 0, 0])
                cube([shim_slot_d*2, shim_slot_w, L+2], center=true);
    }
}

module slice_register() {
    // aneis finos que posicionam cada disco-holder no z correto
    for (k=[0:n_slices-1]) {
        z = -L/2 + slice_thick/2 + k*slice_sep;
        translate([0,0,z]) difference() {
            cylinder(h=2, r=slice_OD/2 + shell_clear, center=true);
            cylinder(h=3, r=slice_OD/2 - 4, center=true);
        }
    }
}

module foot() {
    // base de bancada: trilho sob o cilindro (em -Y), ao longo de todo o eixo z
    r_out = slice_OD/2 + shell_clear + 8;
    translate([0, -(r_out + foot_h/2 - 5), 0])
        cube([foot_w, foot_h, L + 2*plate_t], center=true);
}

// ---------------------------------------------------------------------------
module support() {
    // placas de extremidade
    translate([0,0, plate_z]) end_plate();
    translate([0,0,-plate_z]) rotate([180,0,0]) end_plate();
    tie_rods();
    slice_register();
    shim_former();
    foot();   // base unica de bancada ao longo do eixo
}

support();
