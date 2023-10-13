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
