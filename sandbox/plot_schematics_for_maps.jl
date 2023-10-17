using CubedSphere
using LightningCubedSphere
using GLMakie
using Colors
GLMakie.activate!()


# Figure
fig = Figure(resolution=(2900, 800))

# Zoom
s1 = 1.2
s2 = 0.9
s3 = 1.1
s4 = 1.4

# Axes
ax1 = Axis3(fig[1, 1], aspect = (1, 1, 1), limits = ((-s1, s1), (-s1, s1), (-s1, s1)))
ax2 = Axis3(fig[1, 2], aspect = (1, 1, 1), limits = ((-s2, s2), (-s2, s2), (-s2, s2)))
ax3 = Axis(fig[1, 3], aspect = 1, limits = ((-s3, s3), (-s3, s3)))
ax4 = Axis(fig[1, 4], aspect = 1, limits = ((-s4, s4), (-s4, s4)))
for ax in [ax1, ax2, ax3, ax4]
    hidedecorations!(ax)
    hidespines!(ax)
end


##########
## Cube ##
##########

# Invisible edges, need transparency=false to appear behind
kw = (color=RGBA(0, 0, 0, 1), linewidth=6, transparency=false)
lines!(ax1, [-1, 1], [ 1,  1], [-1, -1]; kw...)
lines!(ax1, [ 1,  1], [-1, 1], [-1, -1]; kw...)
lines!(ax1, [ 1,  1], [ 1,  1], [-1, 1]; kw...)

# Faces
kw = (color=fill(RGBA(0.7, 0.7, 0.7, 0.8), 2, 2), shading=false)
u = [-1 -1;  1 1]
v = [-1  1; -1 1]
w = [ 1  1;  1 1]
surface!(ax1, u, v, -w; kw...) # bottom
surface!(ax1, -w, u', v'; kw...) # side
surface!(ax1, v', -w, u'; kw...) # side
surface!(ax1, u, v, w, color=0*w, shading=false) # top

# Corner dots
kw = (markersize=8, color=:black, transparency=true)
scatter!(ax1,  1,  1,  1; kw...)
scatter!(ax1, -1,  1,  1; kw...)
scatter!(ax1,  1, -1,  1; kw...)
scatter!(ax1, -1, -1,  1; kw...)
scatter!(ax1,  1,  1, -1; kw...)
scatter!(ax1, -1,  1, -1; kw...)
scatter!(ax1,  1, -1, -1; kw...)
scatter!(ax1, -1, -1, -1; kw...)

# Visible edges, need transparency=true to reduce artifacts
kw = (color=RGBA(0, 0, 0, 1), linewidth=6, transparency=true)
lines!(ax1, [-1, 1], [ 1,  1], [ 1,  1]; kw...)
lines!(ax1, [-1, 1], [-1, -1], [ 1,  1]; kw...)
lines!(ax1, [-1, 1], [-1, -1], [-1, -1]; kw...)
lines!(ax1, [ 1,  1], [-1, 1], [ 1,  1]; kw...)
lines!(ax1, [-1, -1], [-1, 1], [ 1,  1]; kw...)
lines!(ax1, [-1, -1], [-1, 1], [-1, -1]; kw...)
lines!(ax1, [-1, -1], [ 1,  1], [-1, 1]; kw...)
lines!(ax1, [ 1,  1], [-1, -1], [-1, 1]; kw...)
lines!(ax1, [-1, -1], [-1, -1], [-1, 1]; kw...)

# North pole dot
scatter!(ax1, 0, 0, 1+0.05; markersize=20, color=:black)


############
## Sphere ##
############

N = 201

# Equator
u = range(0, stop=2π, length=N)
r = 1.01
lines!(ax2, r*cos.(u), r*sin.(u), 0*u, linewidth=6, color=:black)

# conformal_cubed_sphere_mapping
u = range(0, stop=2π, length=N)
v = range(0, stop=π,  length=N)

x = @. cos(u) * sin(v)'
y = @. sin(u) * sin(v)'
z = @. 0u + cos(v)'

surface!(ax2, x , y, z; color=fill(RGBA(0.7, 0.7, 0.7, 0.8), N, N), shading=false)

# Patch
x = range(-1, 1, length=N)
y = range(-1, 1, length=N)
X = x .+ 0*y'
Y = 0*x .+ y'
S = conformal_cubed_sphere_mapping.(X, Y)
sx = getindex.(S, 1)
sy = getindex.(S, 2)
sz = getindex.(S, 3)
surface!(ax2, sx, sy, sz.+1/N, color=0*sz, shading=false)

# North pole dot
scatter!(ax2, 0, 0, 1+0.05, markersize=20, color=:black)


#############
## w plane ##
#############

# Domain
z = s_to_z(sx, sy, sz)
surface!(ax3, real(z), imag(z), 0*sz, shading=false)

# North pole dot
scatter!(ax3, 0, 0, 1, markersize=20, color=:black)


#############
## z plane ##
#############

# Domain
u = [-1 -1;  1 1]
v = [-1  1; -1 1]
surface!(ax4, u, v, 0*u, shading=false)

# North pole dot
scatter!(ax4, 0, 0, 0, markersize=20, color=:black)


# Save
save(joinpath(@__DIR__, "map_schematics.png"), fig)

