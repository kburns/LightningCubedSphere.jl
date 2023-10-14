

function conformal_cubed_sphere_mapping(x, y, W_function)
    X = xᶜ = abs(x)
    Y = yᶜ = abs(y)

    kxy = yᶜ > xᶜ

    xᶜ = 1 - xᶜ
    yᶜ = 1 - yᶜ

    kxy && (xᶜ = 1 - Y)
    kxy && (yᶜ = 1 - X)

    Z = ((xᶜ + im * yᶜ) / 2)^4
    W = W_function(Z)

    im³ = im^(1/3)
    ra = √3 - 1
    cb = -1 + im
    cc = ra * cb / 2

    W = im³ * (W * im)^(1/3)
    W = (W - ra) / (cb + cc * W)
    X, Y = reim(W)

    H = 2 / (1 + X^2 + Y^2)
    X = X * H
    Y = Y * H
    Z = H - 1

    if kxy
        X, Y = Y, X
    end

    y < 0 && (Y = -Y)
    x < 0 && (X = -X)

    # Fix truncation for x = 0 or y = 0.
    x == 0 && (X = 0)
    y == 0 && (Y = 0)

    return X, Y, Z
end

