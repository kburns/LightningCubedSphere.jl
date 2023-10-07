using Printf
using LinearAlgebra
using ForwardDiff
using CubedSphere

"""
    gnomic_projection(cx, cy)

Gnomic projection from northern face of circumscribed cube to northern spherical patch.

Input: (cx, cy) ∈ [-1, 1]²
Output: (sx, sy, sz) ∈ 𝕊²
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

Input: (sx, sy, sz) ∈ 𝕊²
Output: z ∈ ℂ
"""
function stereographic_projection(sx, sy, sz)
    z = @. 2 * (sx + im * sy) / (1 + sz)
    return z
end
SP = stereographic_projection
s_to_z = SP

function inverse_stereographic_projection(z)
    throw("unimplemented")
    sx = sy = sz = nothing
    return (sx, sy, sz)
end
ISP = inverse_stereographic_projection
z_to_s = ISP

"Compose gnomic and stereographic projections."
c_to_z(cx, cy) = SP(GP(cx, cy)...)

"Extend complex vector by 4-fold rotational symmetry."
extend_C4(z) = vcat(z, im*z, -z, -im*z)

"Extend complex vector by degree-4 dihedral symmetry."
extend_D4(z) = extend_C4(vcat(z, reverse(im * conj(z))))

"Sample half-edge (x=1, 0≤y≤1) exponentially close to corner on circumscribed cube."
function sample_edge_tanh(n, w)
    cx = ones(n)
    cy = tanh.(LinRange(0, w, n))
    return (cx, cy)
end

"Sample half-edge (x=1, 0≤y≤1) linearly on circumscribed cube."
function sample_edge_linear(n)
    cx = ones(n)
    cy = LinRange(0, 1, n)
    return (cx, cy)
end

"Samples poles root-exponentially close to a corner point in ℂ."
function sample_poles(n, corner, σ)
    cluster = @. exp(-σ * (sqrt(n) - sqrt(1:n)))
    poles = @. corner * (1 + cluster)
    return poles
end

"runge(z, b) = sum_{j=1}^{n_b} b_j z^(1 + 4(j-1))"
function runge(z, b)
    f = zero(z)
    for (j, bj) in enumerate(b)
        fj = @. bj * z^(1 + 4 * (j-1))
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

"Lightning representation."
lightning(z, p, a, b) = newman(z, p, a) + runge(z, b)

"""Least-squares fit via autodiff: f'(q0) * q = b"""
function fit(f, b, q0)
    # Take Jacobian around initial parameters
    A = ForwardDiff.jacobian(f, q0)
    # Compute column norms
    C = diagm([1 / norm(A[:, j]) for j in axes(A, 2)])
    # Fit
    q = C * ((A * C) \ b)
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
make_plots = true


# Edge samples
c_edge = sample_edge_tanh(ns, w)
c_corner = (1, 1)
z_edge = c_to_z(c_edge...)
z_corner = c_to_z(c_corner...)
@printf "epsilon edge: %.2e (f)\n" minimum(abs.(z_edge .- z_corner))

"""
    reim(z)

Return a tuple with the real and imaginary part of `z`.
"""
reim(z) = real(z), imag(z)

if make_plots
    using GLMakie
    
    # set up the figure
    fig = Figure(resolution = (1600, 1800), fontsize=20)

    kwargs = (aspect = 1, limits = ((-1.7, 1.7), (-1.7, 1.7)))
    ax11 = Axis(fig[1, 1]; kwargs...)
    ax12 = Axis(fig[1, 2]; kwargs...)
    ax13 = Axis(fig[1, 3]; kwargs...)
    ax21 = Axis(fig[2, 1]; kwargs...)
    ax22 = Axis(fig[2, 2]; kwargs...)
    ax23 = Axis(fig[2, 3]; kwargs...)
    ax31 = Axis(fig[3, 1]; kwargs...)
    ax32 = Axis(fig[3, 2]; kwargs...)
    ax33 = Axis(fig[3, 3]; kwargs...)

    for ax in [ax11, ax12, ax13, ax21, ax22, ax23, ax31, ax32, ax33]
        hidedecorations!(ax)
    end
end

if make_plots
    scatter!(ax11, reim.(extend_D4(z_edge)), markersize=8)
end

# Poles
z_poles = sample_poles(na, z_corner, σ)
@printf "epsilon pole: %.2e (f)\n" minimum(abs.(z_poles .- z_corner))
if make_plots
    scatter!(ax11, reim.(extend_C4(z_poles)), markersize=8, color=(:red, 0.8))
end

# Poles
z_poles = sample_poles(na, z_corner, σ)
@printf "epsilon pole: %.2e (f)\n" minimum(abs.(z_poles .- z_corner))
if make_plots
    scatter!(ax11, reim.(extend_C4(z_poles)), markersize=8)
end

# Forward map f: C(stereo) -> C(square)
# Least squares fit: Re(f(z_edge)) = 1
Re_f(q) = real(lightning(z_edge, z_poles, q[1:na], q[na+1:end]))
f_q = fit(Re_f, ones(ns), zeros(na + nb))
println(" (f)")
f_a = f_q[1:na]
f_b = f_q[na+1:end]
forward(z) = lightning(z, z_poles, f_a, f_b)

# Resample and check forward error
c_edge_resample = sample_edge_linear(resample * ns)
z_edge_resample = c_to_z(c_edge_resample...)
y_edge_resample = forward(z_edge_resample)
@printf "resampling error: %.2e (f)\n" maximum(abs.(real(y_edge_resample) .- 1))

# Inverse map g: C(square) -> C(stereo)
# Least squares fit: g(y_edge) = z_edge
y_edge = forward(z_edge)
y_corner = 1 + im
@printf "epsilon edge: %.2e (g)\n" minimum(abs.(y_edge .- y_corner))
y_poles = sample_poles(na, y_corner, σ)
@printf "epsilon pole: %.2e (g)\n" minimum(abs.(y_poles .- y_corner))
Re_g(q) = real(lightning(y_edge, y_poles, q[1:na], q[na+1:end]))
Im_g(q) = imag(lightning(y_edge, y_poles, q[1:na], q[na+1:end]))
ReIm_g(q) = vcat(Re_g(q), Im_g(q))
g_q = fit(ReIm_g, vcat(real(z_edge), imag(z_edge)), zeros(na + nb))
println(" (g)")
g_a = g_q[1:na]
g_b = g_q[na+1:end]
backward(y) = lightning(y, y_poles, g_a, g_b)

# Check backward error
@printf "resampling error: %.2e (g)\n" maximum(abs.(backward(y_edge_resample) - z_edge_resample))

if make_plots
    scatter!(ax12, reim.(extend_D4(y_edge)), markersize=8, color=:green)
    scatter!(ax13, reim.(extend_D4(backward(y_edge))), markersize=8, color=:purple)
end

function plot_gridlines(z_grid, c_grid, y_grid, zz_grid)
    for i in axes(z_grid, 1)
        zᵢ = z_grid[i, :]
        cᵢ = c_grid[i, :]
        yᵢ = y_grid[i, :]
        zzᵢ = zz_grid[i, :]

        lines!(ax21, reim.(cᵢ), color=:blue, linewidth=3)
        lines!(ax22, reim.(zᵢ), color=:black, linewidth=3)
        lines!(ax23, reim.(yᵢ), color=:black, linewidth=3)
        lines!(ax32, reim.(zzᵢ), color=:black, linewidth=3)
        lines!(ax33, reim.(cᵢ), color=:black, linewidth=3)
    end

    for i in axes(z_grid, 2)
        zᵢ = z_grid[:, i]
        cᵢ = c_grid[:, i]
        yᵢ = y_grid[:, i]
        zzᵢ = zz_grid[:, i]

        lines!(ax21, reim.(cᵢ), color=:blue, linewidth=3)
        lines!(ax22, reim.(zᵢ), color=:black, linewidth=3)
        lines!(ax23, reim.(yᵢ), color=:black, linewidth=3)
        lines!(ax32, reim.(zzᵢ), color=:black, linewidth=3)
        lines!(ax33, reim.(cᵢ), color=:black, linewidth=3)
    end
end

function get_various_grids(cx_grid, cy_grid)
    c_grid = @. cx_grid + im * cy_grid
    z_grid = c_to_z(cx_grid, cy_grid)
    y_grid = forward(z_grid)
    zz_grid = backward(c_grid)

    return c_grid, z_grid, y_grid, zz_grid
end

# Plot grids
if make_plots
    cx_grid = range(-1, 1, length = 20)'
    cy_grid = range(-1, 1, length = 20)
    c_grid, z_grid, y_grid, zz_grid = get_various_grids(cx_grid, cy_grid)

    # high_zoom = 1000
    # cx_high = range(-1/high_zoom, 1/high_zoom, length=100)'
    # cy_high = range(-1/high_zoom, 1/high_zoom, length=100)
    # c_grid, z_grid, y_grid, zz_grid = get_various_grids(cx_high, cy_high)

    plot_gridlines(z_grid, c_grid, y_grid, zz_grid)
    
    save("figure.png", fig)
    
    # for i in axes(cz_high, 1)
    #     ci_high = cz_high[i,:]
    #     plot!(10*(9 .+ 3im .+ high_zoom*lightning(ci_high, pInv, aInv, bInv)), color="black", linewidth=0.1)
    # end
    # for i in axes(cz_high, 2)
    #     ci_high = cz_high[:,i]
    #     plot!(10*(9 .+ 3im .+ high_zoom*lightning(ci_high, pInv, aInv, bInv)), color="black", linewidth=0.1)
    # end
    # SPI(X) = SP(X...)
    # cz_high_rancic = SPI.(rancic_s.(cz_high))
    # for i in axes(cz_high_rancic, 1)
    #     ci_high = cz_high_rancic[i,:]
    #     plot!(10*(9 .+ 3im .+ high_zoom*ci_high), color="red", linewidth=0.1)
    # end
    # for i in axes(cz_high_rancic, 2)
    #     ci_high = cz_high_rancic[:,i]
    #     plot!(10*(9 .+ 3im .+ high_zoom*ci_high), color="red", linewidth=0.1)
    # end
    # plot!()
end

if make_plots
    display(fig)
end


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
