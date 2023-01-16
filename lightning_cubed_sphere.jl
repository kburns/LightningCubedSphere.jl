using Printf
using LinearAlgebra
using ForwardDiff
using CubedSphere

"""
    gnomic_projection(cx, cy)

Gnomic projection from northern face of circumscribed cube to northern spherical patch.

Input: (cx, cy) in [-1, 1]^2
Output: (sx, sy, sz) on S^2
"""
function gnomic_projection(cx, cy)
    (maximum(abs.(cx)) > 1 || maximum(abs.(cy)) > 1) && error("both cx and cy must be within [-1, 1]")

    cz = 1
    cr = @. sqrt(cx^2 + cy^2 + cz^2)
    sx = @. cx / cr
    sy = @. cy / cr
    sz = @. cz / cr
    return (sx, sy, sz)
end
GP = gnomic_projection

function inverse_gnomic_projection(sx, sy, sz)
    cx = @. sx / sz
    cy = @. sy / sz
    return (cx, cy)
end
IGP = inverse_gnomic_projection

"""
Stereographic projection from northern spherical patch to C.
Input: (sx, sy, sz) on S^2
Output: z in C
"""
function stereographic_projection(sx, sy, sz)
    z = @. 2 * (sx + im * sy) / (1 + sz)
    return z
end
SP = stereographic_projection

function inverse_stereographic_projection(z)
    throw("unimplemented")
    sx = sy = sz = nothing
    return (sx, sy, sz)
end
ISP = inverse_stereographic_projection

"Compose gnomic and stereographic projections."
project(cx, cy) = SP(GP(cx, cy)...)

"Extend complex vector by 4-fold rotational symmetry."
extend_C4(z) = vcat(z, im*z, -z, -im*z)

"Extend complex vector by degree-4 dihedral symmetry."
extend_D4(z) = extend_C4(vcat(z, reverse(im*conj(z))))

"Sample half-edge exponentially close to corner on circumscribed cube."
function sample_edge_tanh(n, w)
    cx = ones(n)
    cy = tanh.(LinRange(0, w, n))
    return (cx, cy)
end

"Sample half-edge linearly on circumscribed cube."
function sample_edge_linear(n)
    cx = ones(n)
    cy = LinRange(0, 1, n)
    return (cx, cy)
end

"Samples poles root-exponentially close to corners."
function sample_poles(n, corner, σ)
    cluster = @. exp(-σ * (sqrt(n) - sqrt(1:n)))
    poles = @. corner * (1 + cluster)
    return poles
end

"runge(z, b) = sum_{j=1}^{n_b} b_j z^(1+4(j-1))"
function runge(z, b)
    f = zero(z)
    for (j, bj) in enumerate(b)
        fj = @. bj * z^(1 + 4*(j-1))
        f += fj
    end
    return f
end

"newman(z, p, a) = sum_{j=1}^{n_a} sum_{k=0}^{3} (-1)^k i a_j / (z - i^k p_j)"
function newman(z, p, a)
    f = zero(z)
    for (aj, pj) in zip(a, p)
        for k = 0:3
            fjk = @. (-1)^k * im * aj / (z - im^k * pj)
            f += fjk
        end
    end
    return f
end

"reciprocal_log(z, c, s, a) = ..."
function reciprocal_log(z, c, s, a)
    f = zero(z)
    for (aj, sj) in zip(a, s)
        for k = 0:3
            # this step symmetrizes the branch cut
             fjk = @. im^k * cispi(1/4) * aj / (log(-(z/(im^k * c) - 1)) + im*pi - sj)
             f += fjk
             fjk = @. im^k * cispi(1/4) * aj / (log(-(z/(im^k * c) - 1)) - im*pi - sj)
             f += fjk
         end
    end
    return f
end

"Lightning representation."
lightning(z, p, a, b) = newman(z, p, a) + runge(z, b)

"Log-lightning representation."
log_lightning(z, c, s, a, b) = reciprocal_log(z, c, s, a) + runge(z, b)

"""Least-squares fit via autodiff: f'(q0) * q = b"""
function fit(f, b, q0)
    # Take Jacobian around initial parameters
    A = ForwardDiff.jacobian(f, q0)
    # Compute column norms
    C = diagm([1/norm(A[:, j]) for j in axes(A, 2)])
    # Fit
    q = C * ((A*C) \ b)
    # Check error
    e = A * q - b
    @printf "fitting error: %.2e" maximum(abs.(e))
    return q
end

# Parameters
ns = 300       # number of boundary samples per side
na = 80        # number of poles per corner
nb = 10        # number of polynomial terms
w = 17         # boundary sampling tanh width
σ = 4          # pole compaction
resample = 10  # resampling ratio for testing
make_plots = false
include_log = false

# Edge samples
z = project(sample_edge_tanh(ns, w)...)
c = project(1, 1)
@printf "epsilon edge: %.2e (f)\n" minimum(abs.(z .- c))

reim(z) = real(z), imag(z)

if make_plots
    using GLMakie

    fig = Figure()

    kwargs = (aspect = 1, limits = ((-1.7, 1.7), (-1.7, 1.7)))
    ax11 = Axis(fig[1, 1]; kwargs...)
    ax12 = Axis(fig[1, 2]; kwargs...)
    ax13 = Axis(fig[1, 3]; kwargs...)
    ax21 = Axis(fig[2, 1]; kwargs...)
    ax22 = Axis(fig[2, 2]; kwargs...)
    ax23 = Axis(fig[2, 3]; kwargs...)

    for ax in [ax11, ax12, ax13, ax21, ax22, ax23]
        hidedecorations!(ax)
    end

    scatter!(ax11, reim.(extend_D4(z)), markersize=4)
end

# Poles
p = sample_poles(na, c, σ)
@printf "epsilon pole: %.2e (f)\n" minimum(abs.(p .- c))
if make_plots
    scatter!(ax11, reim.(extend_C4(p)), markersize=4)
end

# Least squares fit: Re(f(zE)) = 1
Re_f(q) = real(lightning(z, p, q[1:na], q[na+1:end]))
q = fit(Re_f, ones(ns), zeros(na+nb))
println(" (f)")
a = q[1:na]
b = q[na+1:end]

# Log least squares fit
if include_log
    s = LinRange(-5, 20, na)
    Re_fLog(q) = real(log_lightning(z, c, s, q[1:na], q[na+1:end]))
    qLog = fit(Re_fLog, ones(ns), zeros(na+nb))
    println(" (fLog)")
    aLog = qLog[1:na]
    bLog = qLog[na+1:end]
end

# Resample and check error
z_r = project(sample_edge_linear(resample * ns)...)
f_r = lightning(z_r, p, a, b)
@printf "resampling error: %.2e (f)\n" maximum(abs.(real(f_r) .- 1))
if include_log
    fLog_r = log_lightning(z_r, c, s, aLog, bLog)
    @printf "resampling error: %.2e (fLog)\n" maximum(abs.(real(fLog_r) .- 1))
end

# Inverse least squares fit
zz = project(sample_edge_tanh(ns, w)...)
Z = lightning(zz, p, a, b)
@printf "epsilon edge: %.2e (fInv)\n" minimum(abs.(Z .- (1+im)))
pInv = sample_poles(na, 1+im, σ)
@printf "epsilon pole: %.2e (fInv)\n" minimum(abs.(pInv .- (1+im)))
Re_fInv(q) = real(lightning(Z, pInv, q[1:na], q[na+1:end]))
Im_fInv(q) = imag(lightning(Z, pInv, q[1:na], q[na+1:end]))
ReIm_fInv(q) = vcat(Re_fInv(q), Im_fInv(q))
qInv = fit(ReIm_fInv, vcat(real(zz), imag(zz)), zeros(na+nb))
println(" (fInv)")
aInv = qInv[1:na]
bInv = qInv[na+1:end]

if make_plots
    fInv_r = lightning(f_r, pInv, aInv, bInv)

    scatter!(ax12, reim.(extend_D4(f_r)), markersize=4)
    scatter!(ax13, reim.(extend_D4(fInv_r)), markersize=4)
end

# Plot grids
if make_plots
    cx = range(-1, 1, length=30)'
    cy = range(-1, 1, length=30)
    cz = cx .+ im*cy
    zg = project(cx, cy)

    for i in axes(zg, 1)
        zi = zg[i, :]
        ci = cz[i, :]

        lines!(ax21, reim.(ci), color=(:blue, 0.25))
        lines!(ax21, reim.(zi), color=:black)

        lines!(ax22, reim.(lightning(zi, p, a, b)), color=:black)
        lines!(ax23, reim.(lightning(ci, pInv, aInv, bInv)), color=:black)
    end

    for i in axes(zg, 2)
        zi = zg[:, i]
        ci = cz[:, i]

        lines!(ax21, reim.(ci), color=(:blue, 0.25))
        lines!(ax21, reim.(zi), color=:black)

        lines!(ax22, reim.(lightning(zi, p, a, b)), color=:black)
        lines!(ax23, reim.(lightning(ci, pInv, aInv, bInv)), color=:black)
    end

    display(fig)
    save("figure.png", fig)
end


# Compare to clima's Rancic implementation
rancic_z(X, Y, Z) = (Base.splat(complex) ∘ CubedSphere.conformal_cubed_sphere_inverse_mapping)(X, Y, Z)
rancic_s(z) = CubedSphere.conformal_cubed_sphere_mapping(real(z), imag(z))

# Test on dense linear grid
cx = LinRange(0, 1, ns)'
cy = LinRange(0, 1, ns)
X, Y, Z = GP(cx, cy)
z = SP(X, Y, Z)
z_lightning = lightning(z, p, a, b)
z_rancic = rancic_z.(X, Y, Z)
w_lightning = lightning(z_lightning, pInv, aInv, bInv)
w_rancic = (Base.splat(SP) ∘ rancic_s).(z_rancic)

@printf "lightning roundtrip error: %.2e\n" maximum(abs.(w_lightning - z))
@printf "rancic roundtrip error: %.2e\n" maximum(abs.(w_rancic - z))
@printf "lightning vs rancic forward error: %.2e\n" maximum(abs.(z_lightning - z_rancic))


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
