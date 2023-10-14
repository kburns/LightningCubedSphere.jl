function ReIm_forward(X)
    z = X[1] + im*X[2]
    y = forward(z)
    return vcat(real(y), imag(y))
end

d_f(z) = ForwardDiff.jacobian(ReIm_forward, [real(z), imag(z)])

function ReIm_backward(X)
    y = X[1] + im*X[2]
    z = backward(y)
    return vcat(real(z), imag(z))
end

d_g(y) = ForwardDiff.jacobian(ReIm_backward, [real(y), imag(y)])

cx_dense = LinRange(0, 1, 1000)'
cy_dense = LinRange(0, 1, 1000)
z_dense = c_to_z(cx_dense, cy_dense)
y_dense = forward(z_dense)

J_f = d_f.(z_dense)
J_g = d_g.(y_dense)
JJ_f = norm.(J_f)
JJ_g = norm.(J_g)

