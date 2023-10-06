
using Printf
using Plots
using LinearAlgebra
using ForwardDiff
using CubedSphere

"""
Gnomic projection from northern face of circumscribed cube to northern spherical patch.
Input: (cx, cy) in [-1, 1]^2
Output: (sx, sy, sz) on S^2
"""
function gnomic_projection(cx, cy)
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
    z = @. 2 * (sx + 1im*sy) / (1 + sz)
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
extend_D4(z) = extend_C4(vcat(z, reverse(im*conj(z))))

"Sample half-edge (x=1, 0<=y<=1) exponentially close to corner on circumscribed cube."
function sample_edge_tanh(n, w)
    cx = ones(n)
    cy = tanh.(LinRange(0, w, n))
    return (cx, cy)
end

"Sample half-edge (x=1, 0<=y<=1) linearly on circumscribed cube."
function sample_edge_linear(n)
    cx = ones(n)
    cy = LinRange(0, 1, n)
    return (cx, cy)
end

"Samples poles root-exponentially close to a corner point in C."
function sample_poles(n, corner, σ)
    cluster = @. exp(-σ * (sqrt(n) - sqrt([1:n;])))
    poles = @. corner * (1 + cluster)
    return poles
end

"runge(z, b) = sum_{j=1}^{n_b} b_j z^(1+4(j-1))"
function runge(z, b)
    f = zero(z)
    for (j, bj) in enumerate(b)
        fj = @. bj * z^(1+4*(j-1))
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
    C = diagm([1/norm(A[:,j]) for j in axes(A, 2)])
    # Fit
    q = C * ((A*C) \ b)
    # Check error
    e = A * q - b
    @printf "fitting error: %.2e" maximum(abs.(e))
    return q
end

# Parameters
ns = 300  # number of boundary samples per side
na = 80  # number of poles per corner
nb = 10  # number of polynomial terms
w = 17  # boundary sampling tanh width
σ = 4  # pole compaction
resample = 10  # resampling ratio for testing
make_plots = true

# Edge samples
c_edge = sample_edge_tanh(ns, w)
c_corner = (1, 1)
z_edge = c_to_z(c_edge...)
z_corner = c_to_z(c_corner...)
@printf "epsilon edge: %.2e (f)\n" minimum(abs.(z_edge .- z_corner))
if make_plots
    fig = plot(legend=false, grid=false)
    scatter!(extend_D4(z_edge), aspect_ratio=1, msw=0, ms=1)
end

# Poles
z_poles = sample_poles(na, z_corner, σ)
@printf "epsilon pole: %.2e (f)\n" minimum(abs.(z_poles .- z_corner))
if make_plots
    scatter!(extend_C4(z_poles), aspect_ratio=1, msw=0, ms=1)
end

# Forward map f: C(stereo) -> C(square)
# Least squares fit: Re(f(z_edge)) = 1
Re_f(q) = real(lightning(z_edge, z_poles, q[1:na], q[na+1:end]))
f_q = fit(Re_f, ones(ns), zeros(na+nb))
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
g_q = fit(ReIm_g, vcat(real(z_edge), imag(z_edge)), zeros(na+nb))
println(" (g)")
g_a = g_q[1:na]
g_b = g_q[na+1:end]
backward(y) = lightning(y, y_poles, g_a, g_b)

# Check backward error
@printf "resampling error: %.2e (g)\n" maximum(abs.(backward(y_edge_resample) - z_edge_resample))

if make_plots
    scatter!(3 .+ extend_D4(y_edge), msw=0, ms=1)
    scatter!(6 .+ extend_D4(backward(y_edge)), msw=0, ms=1)
end

# Plot grids
if make_plots
    cx_grid = range(-1, 1, length=20)'
    cy_grid = range(-1, 1, length=20)
    c_grid = cx_grid .+ im*cy_grid
    z_grid = c_to_z(cx_grid, cy_grid)
    y_grid = forward(z_grid)
    zz_grid = backward(c_grid)

    # high_zoom = 1000
    # cx_high = range(-1/high_zoom, 1/high_zoom, length=100)'
    # cy_high = range(-1/high_zoom, 1/high_zoom, length=100)
    # cz_high = cx_high .+ im*cy_high

    for i in axes(z_grid, 1)
        plot!(0 .- 3im .+ c_grid[i,:], color="blue", alpha=0.25)
        plot!(3 .- 3im .+ z_grid[i,:], color="black")
        plot!(6 .- 3im .+ y_grid[i,:], color="black")
        plot!(3 .- 6im .+ zz_grid[i,:], color="black")
        plot!(6 .- 6im .+ c_grid[i,:], color="black")
    end
    for i in axes(z_grid, 2)
        plot!(0 .- 3im .+ c_grid[:,i], color="blue", alpha=0.25)
        plot!(3 .- 3im .+ z_grid[:,i], color="black")
        plot!(6 .- 3im .+ y_grid[:,i], color="black")
        plot!(3 .- 6im .+ zz_grid[:,i], color="black")
        plot!(6 .- 6im .+ c_grid[:,i], color="black")
    end
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
    plot!()
end

# Compare to clima's Rancic implementation
rancic_s_to_y(sx, sy, sz) = (Base.splat(complex) ∘ CubedSphere.conformal_cubed_sphere_inverse_mapping)(sx, sy, sz)
rancic_y_to_s(y) = CubedSphere.conformal_cubed_sphere_mapping(real(y), imag(y))
rancic_forward(z) = rancic_s_to_y.(z_to_s.(z))
#rancic_backward(y) = s_to_z.(rancic_y_to_s.(y))
rancic_backward(y) = (Base.splat(s_to_z) ∘ rancic_y_to_s).(y)

# Test on dense linear grid
cx_dense = LinRange(0, 1, ns)'
cy_dense = LinRange(0, 1, ns)
s_dense = GP(cx_dense, cy_dense)
z_dense = s_to_z(s_dense...)
y_dense = forward(z_dense)
y_dense_rancic = rancic_s_to_y.(s_dense...)
zz_dense = backward(y_dense)
zz_dense_rancic = rancic_backward(y_dense_rancic)
@printf "lightning roundtrip error: %.2e\n" maximum(abs.(zz_dense - z_dense))
@printf "rancic roundtrip error: %.2e\n" maximum(abs.(zz_dense_rancic - z_dense))
@printf "lightning vs rancic forward error: %.2e\n" maximum(abs.(y_dense - y_dense_rancic))

display(fig)

