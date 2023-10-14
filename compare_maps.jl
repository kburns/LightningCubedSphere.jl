using CubedSphere, Printf

include("lightning_test.jl")
include("taylor_coefficients.jl")


# Compare to Clima's Rancic-map implementation

# import W_Rancic, Z_Rancic functions and redefine them using coefficients from taylor_coefficients.jl
import CubedSphere: W_Rancic, Z_Rancic

@info "Use Rancic's original coefficients"

W_Rancic(Z) = sum(A_Rancic[k] * Z^(k-1) for k in length(A_Rancic):-1:1)
Z_Rancic(W) = sum(B_Rancic[k] * W^(k-1) for k in length(B_Rancic):-1:1)

include("compare_lightning_to_rancic_map.jl")


@info "Use MITgcm coefficients as found in the interweb"

A_Rancic = A_MITgcm
B_Rancic = B_MITgcm
W_Rancic(Z) = sum(A_Rancic[k] * Z^(k-1) for k in length(A_Rancic):-1:1)
Z_Rancic(W) = sum(B_Rancic[k] * W^(k-1) for k in length(B_Rancic):-1:1)

include("compare_lightning_to_rancic_map.jl")


@info "Use MITgcm coefficients; A as found in the interweb, B computed"

A_Rancic = A_MITgcm
B_Rancic = B_MITgcm_computed
W_Rancic(Z) = sum(A_Rancic[k] * Z^(k-1) for k in length(A_Rancic):-1:1)
Z_Rancic(W) = sum(B_Rancic[k] * W^(k-1) for k in length(B_Rancic):-1:1)

include("compare_lightning_to_rancic_map.jl")
