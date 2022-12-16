
using Printf
using Plots
using LinearAlgebra
using ForwardDiff

"Gnomic projection from southern face ([-1,1]^2) of unit cube to southern spherical patch."
function gnomic_projection(x, y; R=1)
    z = -1
    r = @. sqrt(x^2 + y^2 + z^2)
    X = @. R * x / r
    Y = @. R * y / r
    Z = @. R * z / r
    return (X, Y, Z)
end
GP = gnomic_projection

"Stereographic projection of southern spherical patch to C."
function stereographic_projection(X, Y, Z; R=1)
    z = @. 2 * R * (X + 1im*Y) / (R - Z)
    return z
end
SP = stereographic_projection

"Extend by 4-fold rotational symmetry."
extend_C4(z) = vcat(z, im*z, -z, -im*z)

"Extend by degree-4 dihedral symmetry."
extend_D4(z) = vcat(extend_C4(z), extend_C4(conj(z)))

"Sample boundary exponentially close to corners on inscribed square."
function sample_boundary_tanh(n, w)
    s = tanh.(collect(LinRange(0, w, n)))
    p1 = ones(n)
    z = SP(GP(p1, s)...)
    return z
end

"Sample boundary linearly on inscribed square."
function sample_boundary_linear(n)
    s = collect(LinRange(0, 1, n))
    p1 = ones(n)
    z = SP(GP(p1, s)...)
    return z
end

"Samples poles root-exponentially close to corners."
function sample_poles(n, σ)
    p = @. exp(-σ * (sqrt(n) - sqrt([1:n;])))
    zC = SP(GP(1, 1)...)
    zP = @. zC + cispi(1/4)*p
    return zP
end

"f_Runge(z, b) = sum_{j=1}^{nb} b_j z^(1+4(j-1))"
function f_Runge(z, b)
    f = zero(z)
    for (j, bj) in enumerate(b)
        fj = @. bj * z^(1+4*(j-1))
        f += fj
    end
    return f
end
#fᵣ

"f_Newman(z, Pz, a) = sum_{j=1}^{na} sum_{k=0}^{3} im a_j (-1)^k / (z - im^k z_j)"
function f_Newman(z, zP, a)
    f = zero(z)
    for (aj, zPj) in zip(a, zP)
        for k = 0:3
            fjk = @. im * aj * (-1)^k / (z - im^k * zPj)
            f += fjk
        end
    end
    return f
end

"f_ReciprocalLog(z, Pz, a) = sum_{j=1}^{nl} sum_{k=0}^{3} a_j e^{i pi/4} i^k / (z - im^k z_j)"
function f_ReciprocalLog(z, sP, a)
    f = zero(z)
    for (aj, sPj) in zip(a, sP)
        for k = 0:3
            fjk = @. aj * exp(im*pi/4) * (im)^k / (log(-(z/(im^k * zC) - 1)) - sPj) # added a negative in front of z so that the branches point away from the square
            f += fjk
        end
    end
    return f
end

"f(z) = f_Newman(z) + f_Runge(z)"
f(z, Pz, a, b) = f_Newman(z, Pz, a) + f_Runge(z, b)

"fl(z) = f_Newman(z) + f_Runge(z)"
fLog(z, Ps, a, b) = f_ReciprocalLog(z, Ps, a) + f_Runge(z, b)

Ps = LinRange(-10,0,na)

# Parameters
ns = 300  # number of boundary samples per side
na = 85  # number of poles per corner
nb = 10  # number of polynomial terms
w = 17  # boundary sampling tanh width
σ = 4  # pole compaction
resample = 10  # resampling ratio for testing

# Compute and plot boundary samples
zE = sample_boundary_tanh(ns, w)
zC = SP(GP(1, 1)...)
@printf "epsilon side: %.2e \n" minimum(abs.(zE .- zC))
plot(extend_D4(zE), aspect_ratio=1, seriestype=:scatter)

# Compute and plot poles
zP = sample_poles(na, σ)
@printf "epsilon pole: %.2e \n" minimum(abs.(zP .- zC))
plot!(extend_C4(zP), aspect_ratio=1, seriestype=:scatter)

# Least squares fit: Re(fE) = 1
Re_fq(q) = real(f(zE, zP, q[1:na], q[na+1:end]))
A = ForwardDiff.jacobian(Re_fq, zeros(na+nb))
B = ones(ns)
C = diagm([1/norm(A[:,j]) for j in axes(A, 2)])  # Normalize columns
Q = C * ((A*C) \ B)
E = A * Q - B
@printf "fitting error: %.2e \n" maximum(abs.(E))

# Reciprocal log least squares fit: Re(fE) = 1
Re_fq_Log(q) = real(fLog(zE, Ps, q[1:na], q[na+1:end]))
ALog = ForwardDiff.jacobian(Re_fq_Log, zeros(na+nb))
BLog = ones(ns)
CLog = diagm([1/norm(ALog[:,j]) for j in axes(ALog, 2)])  # Normalize columns
QLog = CLog * ((ALog*CLog) \ BLog)
ELog = ALog * QLog - BLog
@printf "log fitting error: %.2e \n" maximum(abs.(ELog))

# Extract fitted parameters
a = Q[1:na]
b = Q[na+1:end]

# Extract fitted parameters for reciprocal log
aLog = QLog[1:na]
bLog = QLog[na+1:end]

# Resample and check error
ns_rs = resample * ns
zE_rs = vcat(sample_boundary_linear(ns_rs)...)
fE_rs = f(zE_rs, zP, a, b)
fE_rs_Log = fLog(zE_rs, Ps, aLog, bLog)
EE = real(fE_rs) .- 1
EELog = real(fE_rs_Log) .- 1
@printf "resampled error: %.2e \n" maximum(abs.(EE))
@printf "log resampled error: %.2e \n" maximum(abs.(EELog))

plot!(4 .+ extend_D4(fE_rs), seriestype=:scatter)
plot!(8 .+ extend_D4(fE_rs_Log), seriestype=:scatter)
