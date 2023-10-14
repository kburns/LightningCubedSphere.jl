include("LightningMap.jl")

using .LightningMap

# Parameters
na = 80        # number of poles per corner
nb = 10        # number of polynomial terms
ns = 300       # number of boundary samples per side
w = 17         # boundary sampling tanh width
σ = 4          # pole compaction (Gopal & Trefethen (SINUM 2019))
resample = 4   # resampling ratio for testing

forward, backward = LightningMap.compute_lightning_maps(ns, na, nb, w, σ, resample);
