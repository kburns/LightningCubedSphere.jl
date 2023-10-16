using CubedSphere
using LightningCubedSphere
using GLMakie
using Colors
GLMakie.activate!()

fig = Figure(resolution=(2900, 800))

ax1 = Axis3(fig[1, 1], aspect = (1, 1, 1), limits = ((-1.3, 1.3), (-1.3, 1.3), (-1.3, 1.3)))
ax2 = Axis3(fig[1, 2], aspect = (1, 1, 1), limits = ((-1.05, 1.05), (-1.05, 1.05), (-1.05, 1.05)))
ax3 = Axis(fig[1, 3], aspect = 1, limits = ((-1, 1), (-1, 1)))
ax4 = Axis(fig[1, 4], aspect = 1, limits = ((-1.5, 1.5), (-1.5, 1.5)))

for ax in [ax1, ax2, ax3, ax4]
    hidedecorations!(ax)
    hidespines!(ax)
end

# the conformal map domain

N = 256

x = range(-1, 1, length=N)
y = range(-1, 1, length=N)

X = zeros(length(x), length(y))
Y = zeros(length(x), length(y))
Z = zeros(length(x), length(y))

for (j, y′) in enumerate(y), (i, x′) in enumerate(x)
    X[i, j], Y[i, j], Z[i, j] = conformal_cubed_sphere_mapping(x′, y′)
end

surface!(ax2, X, Y, Z, color=0*Z, alpha=0.8, linewidth=4)

surface!(ax3, X, Y, 0*Z, shading = false)


# Equator on sphere

N = 256

x = range(-1, 1, length=N)
y = range(-1, 1, length=N)
Zc = zeros(length(x), length(y))
[Zc[i, j] = x[i]^2 + y[j]^2 for i in 1:length(x), j in 1:length(y)]

contour!(ax2, x, y, Zc, levels=[1.00001], color=RGBA(0.1, 0.1, 0.1, 1), linewidth=6)


# Sphere

n = 101
u = range(0, stop=2π, length=n)
v = range(0, stop=π,  length=n)

x = zeros(n, n)
y = zeros(n, n)
z = zeros(n, n)

R = 1

for i in 1:n
    for j in 1:n
        x[i, j] = R * cos(u[i]) * sin(v[j])
        y[i, j] = R * sin(u[i]) * sin(v[j])
        z[i, j] = R * cos(v[j])
    end
end

surface!(ax2, x , y, z, color = fill(RGBA(0.7, 0.7, 0.7, 0.8), n, n), shading=false)


# Domain on cube

n = 101
u = range(-1, 1, length=n)
v = range(-1, 1,  length=n)

X = zeros(length(u), length(v))
Y = zeros(length(u), length(v))
Z = zeros(length(u), length(v))

for i in 1:n, j in 1:n
    X[i, j], Y[i, j], Z[i, j] = (u[i], v[j], 1)
end

surface!(ax1, X, Y, Z, color=0*Z, alpha=0.8, linewidth=4, shading = true)


# Domain on panel (d)

surface!(ax4, X, Y, 0*X .+ 1, shading = false)


# Cube's edges

lines!(ax1, [-1, 1], [1, 1], [1, 1], color=RGBA(0.1, 0.1, 0.1, 1), linewidth=6)
lines!(ax1, [-1, 1], [-1, -1], [1, 1], color=RGBA(0.1, 0.1, 0.1, 1), linewidth=6)
lines!(ax1, [-1, 1], [1, 1], [-1, -1], color=RGBA(0.1, 0.1, 0.1, 1), linewidth=6)
lines!(ax1, [-1, 1], [-1, -1], [-1, -1], color=RGBA(0.1, 0.1, 0.1, 1), linewidth=6)

lines!(ax1, [1,   1], [-1, 1], [1, 1], color=RGBA(0.1, 0.1, 0.1, 1), linewidth=6)
lines!(ax1, [-1, -1], [-1, 1], [1, 1], color=RGBA(0.1, 0.1, 0.1, 1), linewidth=6)
lines!(ax1, [1,   1], [-1, 1], [-1, -1], color=RGBA(0.1, 0.1, 0.1, 1), linewidth=6)
lines!(ax1, [-1, -1], [-1, 1], [-1, -1], color=RGBA(0.1, 0.1, 0.1, 1), linewidth=6)

lines!(ax1, [1,   1], [ 1,  1], [-1, 1], color=RGBA(0.1, 0.1, 0.1, 1), linewidth=6)
lines!(ax1, [-1, -1], [ 1,  1], [-1, 1], color=RGBA(0.1, 0.1, 0.1, 1), linewidth=6)
lines!(ax1, [1,   1], [-1, -1], [-1, 1], color=RGBA(0.1, 0.1, 0.1, 1), linewidth=6)
lines!(ax1, [-1, -1], [-1, -1], [-1, 1], color=RGBA(0.1, 0.1, 0.1, 1), linewidth=6)


# Cube's faces

shading = false
color = fill(RGBA(0.7, 0.7, 0.7, 0.8), n, n)

x = zeros(n, n)
y = zeros(n, n)
z = zeros(n, n)

for i in 1:n, j in 1:n
    x[i, j] = u[i]
    y[i, j] = v[j]
    z[i, j] = -1
end
surface!(ax1, x, y, z; color, shading)

for i in 1:n, j in 1:n
    x[i, j] = -1
    y[i, j] = u[i]
    z[i, j] = v[j]
end
surface!(ax1, x, y, z; color, shading)

for i in 1:n, j in 1:n
    x[i, j] = u[i]
    y[i, j] = -1
    z[i, j] = v[j]
end
surface!(ax1, x, y, z; color, shading)


# North Pole dots

scatter!(ax1, 0, 0, 1, markersize=20, color=:black)
scatter!(ax2, 0, 0, 1, markersize=20, color=:black)
scatter!(ax3, 0, 0, 1, markersize=20, color=:black)
scatter!(ax4, 0, 0, 1, markersize=20, color=:black)

colgap!(fig.layout, 1, Relative(0))

fig

save("maps_domains.png", fig)
