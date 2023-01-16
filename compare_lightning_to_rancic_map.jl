rancic_z(X, Y, Z) = (Base.splat(complex) ∘ conformal_cubed_sphere_inverse_mapping)(X, Y, Z)
rancic_s(z) = conformal_cubed_sphere_mapping(real(z), imag(z))

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
@printf "rancic-method roundtrip error: %.2e\n" maximum(abs.(w_rancic - z))
@printf "lightning vs rancic-method forward error: %.2e\n" maximum(abs.(z_lightning - z_rancic))
