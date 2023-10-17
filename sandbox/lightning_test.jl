
using Printf
using LightningCubedSphere


# Parameters
na = 80        # number of poles per corner
nb = 20        # number of polynomial terms
ns = 300       # number of boundary samples per side
σ = 4          # pole compaction (Gopal & Trefethen (SINUM 2019))
resample = 4   # resampling ratio for testing

forward, backward = compute_lightning_maps(ns, na, nb, σ, resample);
println()

# Test vs Rancic on dense grid
#cx, cy = sample_edge_rootexp(ns, σ/sqrt(ns/na))
cx, cy = sample_edge_linear(ns)
cx = cy'
z = c_to_z(cx, cy)
w = cx .+ im*cy

@printf "Lightning forward+backward error: %.2e\n" maximum(abs.(backward(forward(z)) - z))
@printf "Lightning backward+forward error: %.2e\n" maximum(abs.(forward(backward(w)) - w))
@printf "Rancic forward+backward error: %.2e\n" maximum(abs.(Rancic_backward.(Rancic_forward.(z)) - z))
@printf "Rancic backward+forward error: %.2e\n" maximum(abs.(Rancic_forward.(Rancic_backward.(w)) - w))
@printf "Lightning vs Rancic forward error: %.2e\n" maximum(abs.(forward(z) - Rancic_forward.(z)))
@printf "Lightning vs Rancic backward error: %.2e\n" maximum(abs.(backward(w) - Rancic_backward.(w)))

