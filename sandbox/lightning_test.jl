
using Printf
using LightningCubedSphere


# Parameters
na = 80        # number of poles per corner
nb = 10        # number of polynomial terms
ns = 300       # number of boundary samples per side
σ = 4          # pole compaction (Gopal & Trefethen (SINUM 2019))
resample = 4   # resampling ratio for testing

forward, backward = compute_lightning_maps(ns, na, nb, σ, resample);
println()

# Test on dense linear grid
cx_dense = LinRange(0, 1, ns)'
cy_dense = LinRange(0, 1, ns)
s_dense = c_to_s(cx_dense, cy_dense)
z_dense = s_to_z(s_dense...)
w_dense = forward(z_dense)
w_Rancic = Rancic_forward.(z_dense)
zz_dense = backward(w_dense)
zz_Rancic = Rancic_backward.(w_Rancic)

@printf "lightning roundtrip error: %.2e\n" maximum(abs.(zz_dense - z_dense))
@printf "Rancic roundtrip error: %.2e\n" maximum(abs.(zz_Rancic - z_dense))
@printf "lightning vs Rancic forward error: %.2e\n" maximum(abs.(w_dense - w_Rancic))

w_dense = cx_dense .+ im*cy_dense
zz_dense = backward(w_dense)
zz_Rancic = Rancic_backward.(w_dense)
@printf "lightning vs Rancic backward error: %.2e\n" maximum(abs.(zz_dense - zz_Rancic))

