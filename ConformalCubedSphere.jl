
using Printf
using Plots
using LinearAlgebra
using ForwardDiff

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

"Sample edge exponentially close to corners on circumscribed cube and project to C."
function sample_edge_tanh(n, w)
    cx = ones(n)
    cy = tanh.(LinRange(0, w, n))
    z = project(cx, cy)
    return z
end

"Sample edge linearly on circumscribed cube and project to C."
function sample_edge_linear(n)
    cx = ones(n)
    cy = LinRange(0, 1, n)
    z = project(cx, cy)
    return z
end

"Samples poles root-exponentially close to corners."
function sample_poles(n, σ)
    cluster = @. exp(-σ * (sqrt(n) - sqrt([1:n;])))
    corner = project(1, 1)
    poles = @. corner + cispi(1/4) * cluster
    return poles
end

"f_Runge(z, b) = sum_{j=1}^{n_b} b_j z^(1+4(j-1))"
function f_Runge(z, b)
    f = zero(z)
    for (j, bj) in enumerate(b)
        fj = @. bj * z^(1+4*(j-1))
        f += fj
    end
    return f
end

"f_Newman(z, p, a) = sum_{j=1}^{n_a} sum_{k=0}^{3} (-1)^k i a_j / (z - i^k p_j)"
function f_Newman(z, p, a)
    f = zero(z)
    for (aj, pj) in zip(a, p)
        for k = 0:3
            fjk = @. (-1)^k * im * aj / (z - im^k * pj)
            f += fjk
        end
    end
    return f
end

"f_ReciprocalLog(z, c, s, a) = ..."
function f_ReciprocalLog(z, c, s, a)
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
f(z, p, a, b) = f_Newman(z, p, a) + f_Runge(z, b)

"Log-lightning representation."
fLog(z, c, s, a, b) = f_ReciprocalLog(z, c, s, a) + f_Runge(z, b)

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
na = 85  # number of poles per corner
nb = 10  # number of polynomial terms
w = 17  # boundary sampling tanh width
σ = 4  # pole compaction
resample = 10  # resampling ratio for testing

# Compute and plot edge samples
z = sample_edge_tanh(ns, w)
c = project(1, 1)
@printf "epsilon edge: %.2e \n" minimum(abs.(z .- c))
plot(legend=false, grid=false)
scatter!(extend_D4(z), aspect_ratio=1, msw=0, ms=2)

# Compute and plot poles
p = sample_poles(na, σ)
@printf "epsilon pole: %.2e \n" minimum(abs.(p .- c))
scatter!(extend_C4(p), aspect_ratio=1, msw=0, ms=2)

# Least squares fit: Re(f(zE)) = 1
Re_f(q) = real(f(z, p, q[1:na], q[na+1:end]))
q = fit(Re_f, ones(ns), zeros(na+nb))
println(" (f)")
a = q[1:na]
b = q[na+1:end]

s = LinRange(-5,20,na)
Re_fLog(q) = real(fLog(z, c, s, q[1:na], q[na+1:end]))
qLog = fit(Re_fLog, ones(ns), zeros(na+nb))
println(" (fLog)")
aLog = qLog[1:na]
bLog = qLog[na+1:end]

# Resample and check error
z_r = vcat(sample_edge_linear(resample * ns)...)
f_r = f(z_r, p, a, b)
fLog_r = fLog(z_r, c, s, aLog, bLog)
@printf "resampling error: %.2e (f)\n" maximum(abs.(real(f_r) .- 1))
@printf "resampling error: %.2e (fLog)\n" maximum(abs.(real(fLog_r) .- 1))
scatter!(3 .+ extend_D4(f_r), msw=0, ms=2)
scatter!(6 .+ extend_D4(fLog_r), msw=0, ms=2)

"Inverse lightning representation."
fInv(z, p, a, b) = f_Newman(z, p./abs(c)*abs(1+im), a) + f_Runge(z, b)

# Get real and imaginary parts of the Jacobian since ForwardDiff doesn't work for complex functions
Z = f(z, p, a, b)
Re_inv(q) = real(fInv(Z, p, q[1:na], q[na+1:end]))
Im_inv(q) = imag(fInv(Z, p, q[1:na], q[na+1:end]))

# Take Jacobian around initial parameters
A_inv = ForwardDiff.jacobian(Re_inv, q) + im*ForwardDiff.jacobian(Im_inv, q) 
# Compute column norms
C_inv = diagm([1/norm(A_inv[:,j]) for j in axes(A_inv, 2)])
# Fit
RHS = z
q_inv = C_inv * ((A_inv*C_inv) \ RHS)
# Check error
e = A_inv * q_inv - RHS
@printf "fitting error: %.2e" maximum(abs.(e))
println(" (inv_f)")
a_inv = q_inv[1:na]
b_inv = q_inv[na+1:end]

# Plot grids
cx = range(-1, 1, length=20)'
cy = range(-1, 1, length=20)
cz = cx .+ im*cy
zg = project(cx, cy)
for i in axes(zg, 1)
    zi = zg[i,:]
    ci = cz[i,:]
    plot!(3im .+ cz[i,:], color="blue", alpha=0.25)
    plot!(3im .+ zi, color="black")
    plot!(3 .+ 3im .+ f(zi, p, a, b), color="black")
    di = f(zi, p, a, b)
    plot!(6 .+ 3im .+ fLog(zi, c, s, aLog, bLog), color="black")
    plot!(9 .+ 3im .+ fInv(di, p, a_inv, b_inv), color="black")
end
for i in axes(zg, 2)
    zi = zg[:,i]
    ci = cz[:,i]
    plot!(3im .+ cz[:,i], color="blue", alpha=0.25)
    plot!(3im .+ zi, color="black")
    plot!(3 .+ 3im .+ f(zi, p, a, b), color="black")
    di = f(zi, p, a, b)
    plot!(6 .+ 3im .+ fLog(zi, c, s, aLog, bLog), color="black")
    plot!(9 .+ 3im .+ fInv(di, p, a_inv, b_inv), color="black")
end
plot!()