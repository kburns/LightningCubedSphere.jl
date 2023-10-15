using TaylorSeries

# from Rancic et al. 1996 paper
A_Rancic = [
    +0.00000000000000,
    +1.47713062600964,
    -0.38183510510174,
    -0.05573058001191,
    -0.00895883606818,
    -0.00791315785221,
    -0.00486625437708,
    -0.00329251751279,
    -0.00235481488325,
    -0.00175870527475,
    -0.00135681133278,
    -0.00107459847699,
    -0.00086944475948,
    -0.00071607115121,
    -0.00059867100093,
    -0.00050699063239,
    -0.00043415191279,
    -0.00037541003286,
    -0.00032741060100,
    -0.00028773091482,
    -0.00025458777519,
    -0.00022664642371,
    -0.00020289261022,
    -0.00018254510830,
    -0.00016499474461,
    -0.00014976117168,
    -0.00013646173946,
    -0.00012478875823,
    -0.00011449267279,
    -0.00010536946150,
    -0.00009725109376
]

# from http://wwwcvs.mitgcm.org/viewvc/MITgcm/MITgcm_contrib/high_res_cube/matlab-grid-generator/map_xy2xyz.m?view=markup#l73
A_MITgcm = [
    +0.00000000000000,
    +1.47713057321600,
    -0.38183513110512,
    -0.05573055466344,
    -0.01895884801823,
    -0.00791314396396,
    -0.00486626515498,
    -0.00329250387158,
    -0.00235482619663,
    -0.00175869000970,
    -0.00135682443774,
    -0.00107458043205,
    -0.00086946107050,
    -0.00071604933286,
    -0.00059869243613,
    -0.00050696402446,
    -0.00043418115349,
    -0.00037537743098,
    -0.00032745130951,
    -0.00028769063795,
    -0.00025464473946,
    -0.00022659577923,
    -0.00020297175587,
    -0.00018247947703,
    -0.00016510295548,
    -0.00014967258633,
    -0.00013660647356,
    -0.00012466390509,
    -0.00011468147908,
    -0.00010518717478,
    -0.00009749136078
]

# from http://wwwcvs.mitgcm.org/viewvc/MITgcm/MITgcm_contrib/high_res_cube/matlab-grid-generator/map_xy2xyz.m?view=markup#l73
B_MITgcm = [
    +0.00000000000000,
    +0.67698822171341,
    +0.11847295533659,
    +0.05317179075349,
    +0.02965811274764,
    +0.01912447871071,
    +0.01342566129383,
    +0.00998873721022,
    +0.00774869352561,
    +0.00620347278164,
    +0.00509011141874,
    +0.00425981415542,
    +0.00362309163280,
    +0.00312341651697,
    +0.00272361113245,
    +0.00239838233411,
    +0.00213002038153,
    +0.00190581436893,
    +0.00171644267546,
    +0.00155493871562,
    +0.00141600812949,
    +0.00129556691848,
    +0.00119042232809,
    +0.00109804804853,
    +0.00101642312253,
    +0.00094391466713,
    +0.00087919127990,
    +0.00082115825576,
    +0.00076890854394,
    +0.00072168520663,
    +0.00067885239089
]

# Correct inverse of Rancic coeffs
A_series = Taylor1(A_Rancic)
B_series = inverse(A_series)
B_Rancic = B_series.coeffs

# Correct inverse of MITgcm coeffs
A_series = Taylor1(A_MITgcm)
B_series = inverse(A_series)
B_MITgcm_computed = B_series.coeffs

# Rancic and MITgcm maps
W_Rancic(Z) = sum(A_Rancic[k] * Z^(k-1) for k in length(A_Rancic):-1:1)
W_MITcgm(Z) = sum(A_MITcgm[k] * Z^(k-1) for k in length(A_MITcgm):-1:1)
Z_Rancic(W) = sum(B_Rancic[k] * W^(k-1) for k in length(B_Rancic):-1:1)
Z_MITcgm(W) = sum(B_MITcgm[k] * W^(k-1) for k in length(B_MITcgm):-1:1)
Rancic_y_to_s(y) = conformal_cubed_sphere_mapping(real(y), imag(y); W_map=W_Rancic)
MITgcm_y_to_s(y) = conformal_cubed_sphere_mapping(real(y), imag(y); W_map=W_MITcgm)
Rancic_s_to_y(s) = conformal_cubed_sphere_inverse_mapping(s...; Z_map=Z_Rancic)
MITgcm_s_to_y(s) = conformal_cubed_sphere_inverse_mapping(s...; Z_map=Z_MITcgm)
Rancic_backward(y) = s_to_z(Rancic_y_to_s(y)...)
MITgcm_backward(y) = s_to_z(MITgcm_y_to_s(y)...)
Rancic_forward(z) = Rancic_s_to_y(z_to_s(z)...)
MITgcm_forward(z) = MITgcm_s_to_y(z_to_s(z)...)

