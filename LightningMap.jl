module LightningMap

using ExportAll
using Printf
using LinearAlgebra
using ForwardDiff
using CubedSphere
using CairoMakie


"Gnomic projection from (cx, cy) ∈ [-1, 1]² to northern patch on unit sphere."
function c_to_s(cx, cy)
    cz = 1
    cr = @. sqrt(cx^2 + cy^2 + cz^2)
    sx = @. cx / cr
    sy = @. cy / cr
    sz = @. cz / cr
    return (sx, sy, sz)
end


"Inverse gnomic projection from northern patch on unit sphere to (cx, cy) ∈ [-1, 1]²."
function s_to_c(sx, sy, sz)
    cx = @. sx / sz
    cy = @. sy / sz
    return (cx, cy)
end


"Stereographic projection from unit sphere to northern tangent plane (complex)."
function s_to_z(sx, sy, sz)
    zx = @. 2 * sx / (1 + sz)
    zy = @. 2 * sy / (1 + sz)
    z = @. zx + im * zy
    return z
end


"Inverse stereographic projection from northern tangent plane (complex) to unit sphere."
function z_to_s(z)
    zx = real(z)
    zy = imag(z)
    R2 = @. zx^2 + zy^2
    sz = @. (4 - R2) / (4 + R2)
    sx = @. zx * (1 + sz) / 2
    sy = @. zy * (1 + sz) / 2
    return (sx, sy, sz)
end


"Compose gnomic and stereographic projections."
c_to_z(cx, cy) = s_to_z(c_to_s(cx, cy)...)
z_to_c(z) = s_to_c(z_to_s(z)...)


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


"Sample half-edge (x=1, 0≤y≤1) root-exponentially."
function sample_edge_rootexp(n, σ)
    # Eq 3.2 from Gopal & Trefethen (SINUM 2019)
    cluster = @. exp(-σ * (sqrt(n) - sqrt(1:n)))
    cx = ones(n)
    cy = 1 .- cluster
    return (cx, cy)
end


"Samples poles root-exponentially close to a corner point in ℂ."
function sample_poles(n, corner, σ)
    # Eq 3.2 from Gopal & Trefethen (SINUM 2019)
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


"""Least-squares fit via autodiff: f(q) = b"""
function fit(f, b, q0; iters=1)
    q = q0
    for i in 1:iters
        # Take Jacobian around current parameters
        A = ForwardDiff.jacobian(f, q)
        # Compute column norms
        C = diagm([1 / norm(A[:, j]) for j in axes(A, 2)])
        # Fit
        dq = C * ((A * C) \ (b - f(q)))
        q += dq
    end
    # Check error
    e = f(q) - b
    @printf "fitting error: %.2e" maximum(abs.(e))
    return q
end


"Return a tuple with the real and imaginary part of `z`."
reim(z) = real(z), imag(z)


"""
    compute_lightning_maps(ns, na, nb, w, σ, resample; make_plots=false)

Compute forward and backward conformal maps using the lightning representation.

Input:
- `ns`: number of boundary samples per side
- `na`: number of poles per corner
- `nb`: number of polynomial terms
- `σ`: pole compaction
- `resample`: resampling ratio for testing
- `make_plots`: whether to make plots
"""
function compute_lightning_maps(ns, na, nb, σ, resample; make_plots=false)

    # Edge samples
    c_edge = sample_edge_rootexp(ns, σ/sqrt(ns/na))
    c_corner = (1, 1)
    z_edge = c_to_z(c_edge...)
    z_corner = c_to_z(c_corner...)
    @printf "epsilon edge: %.2e (f)\n" minimum(abs.(z_edge .- z_corner))

    # Poles
    z_poles = sample_poles(na, z_corner, σ)
    @printf "epsilon pole: %.2e (f)\n" minimum(abs.(z_poles .- z_corner))

    # Forward map f: C(stereo) -> C(square)
    # Least squares fit: Re(f(z_edge)) = 1
    Re_f(q) = real(lightning(z_edge, z_poles, q[1:na], q[na+1:end]))
    f_q = fit(Re_f, ones(ns), zeros(na + nb), iters=2)
    println(" (f)")
    f_a = f_q[1:na]
    f_b = f_q[na+1:end]
    forward(z) = lightning(z, z_poles, f_a, f_b)

    # Resample and check forward error
    #c_edge_resample = sample_edge_linear(resample * ns)
    c_edge_resample = sample_edge_rootexp(resample*ns, σ/sqrt(resample*ns/na))
    z_edge_resample = c_to_z(c_edge_resample...)
    y_edge_resample = forward(z_edge_resample)
    @printf "resampling error: %.2e (f)\n" maximum(abs.(real(y_edge_resample) .- 1))
    @printf "corner error: %.2e (f)\n" abs.(forward(z_corner) - (1+im))

    # Inverse map g: C(square) -> C(stereo)
    # Least squares fit: g(f(z_edge)) = z_edge
    y_edge = forward(z_edge)
    y_corner = 1 + im
    @printf "epsilon edge: %.2e (g0)\n" minimum(abs.(y_edge .- y_corner))
    y_poles = sample_poles(na, y_corner, σ)
    @printf "epsilon pole: %.2e (g0)\n" minimum(abs.(y_poles .- y_corner))
    Re_g(q) = real(lightning(y_edge, y_poles, q[1:na], q[na+1:end]))
    Im_g(q) = imag(lightning(y_edge, y_poles, q[1:na], q[na+1:end]))
    ReIm_g(q) = vcat(Re_g(q), Im_g(q))
    g_q = fit(ReIm_g, vcat(real(z_edge), imag(z_edge)), zeros(na + nb))
    println(" (g0)")
    g_a = g_q[1:na]
    g_b = g_q[na+1:end]
    backward(y) = lightning(y, y_poles, g_a, g_b)

    # Check backward error
    @printf "resampling error: %.2e (g0)\n" maximum(abs.(backward(y_edge_resample) - z_edge_resample))
    @printf "corner error: %.2e (g0)\n" abs.(backward(1+im) - z_corner)

    # Inverse map g: C(square) -> C(stereo)
    # Refine via Newton's method: z_to_c(g(y_edge))[1] = 1
    y_edge = sample_edge_rootexp(ns, σ/sqrt(ns/na))
    y_edge = y_edge[1] + im*y_edge[2]
    y_corner = 1 + im
    @printf "epsilon edge: %.2e (g)\n" minimum(abs.(y_edge .- y_corner))
    y_poles = sample_poles(na, y_corner, σ)
    @printf "epsilon pole: %.2e (g)\n" minimum(abs.(y_poles .- y_corner))
    y_to_cx(q) = z_to_c(lightning(y_edge, y_poles, q[1:na], q[na+1:end]))[1]
    g_q = fit(y_to_cx, ones(ns), g_q, iters=2)
    println(" (g)")
    g_a = g_q[1:na]
    g_b = g_q[na+1:end]
    backward(y) = lightning(y, y_poles, g_a, g_b)

    # Check backward error
    @printf "resampling error: %.2e (g)\n" maximum(abs.(backward(y_edge_resample) - z_edge_resample))
    @printf "corner error: %.2e (g)\n" abs.(backward(1+im) - z_corner)

    if make_plots
        # set up the figure
        fig = Figure(resolution = (1600, 600), fontsize=20)

        kwargs = (aspect = 1, limits = ((-1.7, 1.7), (-1.7, 1.7)))
        ax11 = Axis(fig[1, 1]; kwargs...)
        ax12 = Axis(fig[1, 2]; kwargs...)
        ax13 = Axis(fig[1, 3]; kwargs...)

        for ax in [ax11, ax12, ax13]
            hidedecorations!(ax)
        end
        scatter!(ax11, reim.(extend_C4(z_poles)), markersize=8, color=(:red, 0.8))
        scatter!(ax11, reim.(extend_C4(z_poles)), markersize=8)
        scatter!(ax11, reim.(extend_D4(z_edge)), markersize=8)

        scatter!(ax12, reim.(extend_D4(y_edge)), markersize=8, color=:green)
        scatter!(ax13, reim.(extend_D4(backward(y_edge))), markersize=8, color=:purple)

        save("lightning_map.pdf", fig)
    end

    return forward, backward
end

@exportAll()

end # module
