using LightningCubedSphere
using CairoMakie
using CubedSphere


R_earth_km = 6371

# Lightning maps
ns = 300       # number of boundary samples per side
na = 80        # number of poles per corner
nb = 10        # number of polynomial terms
σ = 4          # pole compaction
resample = 10  # resampling ratio for testing
forward, backward = compute_lightning_maps(ns, na, nb, σ, resample; make_plots=true);


"plot the grid lines for a complex-valued `z_grid`"
function plot_gridlines!(z_grid; color=:black, linewidth=2, ax=nothing)
    for i in axes(z_grid, 1)
        zᵢ = z_grid[i, :]
        lines!(ax, reim.(zᵢ); color, linewidth)
    end
    for i in axes(z_grid, 2)
        zᵢ = z_grid[:, i]
        lines!(ax, reim.(zᵢ); color, linewidth)
    end
    return ax
end

function plot_map(w_grid, z_grid; filename=nothing, fig=nothing, ax=nothing, color=:black)
    if isnothing(fig)
        fig = Figure(resolution = (800, 400), fontsize=20)
    end
    if isnothing(ax)
        ax1 = Axis(fig[1, 1:2]; title="llll", aspect=DataAspect())
        # ax2 = Axis(fig[1, 2]; aspect=DataAspect())
        ax = (ax1, ax2)
    end
    plot_gridlines!(w_grid; color, ax=ax[1])
    plot_gridlines!(z_grid; color, ax=ax[2])
    if !isnothing(filename)
        save(filename, fig)
    end
    return fig, ax
end

# Plot maps over full square
wx_grid = range(-1, 1, length=20)'
wy_grid = range(-1, 1, length=20)
w_grid = @. wx_grid + im * wy_grid

# fig, ax = plot_map(w_grid, backward(w_grid); color=:green)
# fig, ax = plot_map(w_grid, mitgcm_backward(w_grid); fig=fig, ax=ax, color=(:red, 0.5), filename="full_map.pdf")


# Plot quadrants
R = R_earth_km
J = 0.8133665002897938

function plot_z_quadrants!(L, N, backward, ax; box=false)

    # Build w grid
    wx_grid = range(-L, L, length=N)
    wy_grid = range(-L, L, length=N)
    w_grid = @. wx_grid' + im * wy_grid
    wx_half_grid = wx_grid[floor(Int,N/2+1):end]
    wy_half_grid = wy_grid[floor(Int,N/2+1):end]
    w_half_grid = @. wx_half_grid' + im * wy_half_grid

    # Evaluate z grids
    z_grid = backward.(w_grid)
    z_half_grid = backward.(w_half_grid)

    # # Rescale grids
    # Z_grid = R * z_grid
    # Z_half_grid = R * z_half_grid

    # Plot grids
    plot_gridlines!(         z_grid; ax=ax, color=(:black,0.4))
    plot_gridlines!(    z_half_grid; ax=ax, color=(:green,0.5))
    plot_gridlines!(-   z_half_grid; ax=ax, color=(:red,0.5))
    plot_gridlines!( im*z_half_grid; ax=ax, color=(:blue,0.5))
    plot_gridlines!(-im*z_half_grid; ax=ax, color=(:purple,0.5))

    # # Add scale bar
    # dx = 2 * L / (N-1) * R * J
    # y0 = minimum(imag(Z_grid)) - dx
    # lines!(ax[2], [-dx/2, dx/2], [y0, y0]; color=:black, linewidth=2)
    # if dx < 1
    #     label = "$(round(Int,dx*1000)) m"
    # else
    #     label = "$(round(Int,dx)) km"
    # end
    # text!(ax[2], 0, y0-dx/2; text=label, align=(:center, :top), color=:black, fontsize=20)

    # Set limits
    # L * R * J
    # limit = 1.2 * R * J * L
    #xlims!(ax, (-limit, limit))
    #ylims!(ax, (-limit, limit))

    # Remove axes
    hidexdecorations!(ax)
    hideydecorations!(ax)

    if box
        xlims!(ax, (minimum(real(z_grid)), maximum(real(z_grid))))
        ylims!(ax, (minimum(imag(z_grid)), maximum(imag(z_grid))))
    else
        hidespines!(ax)
    end

end

fig = Figure(resolution = (1200, 800), fontsize=20)
ax = Axis(fig[1:9, 1:9]; aspect=DataAspect())
plot_z_quadrants!(1, 100, Rancic_backward, ax)
s = 0.02
lines!(ax, [s, s, -s, -s, s], [-s, s, s, -s, -s]; color=:black, linewidth=4)

dx = 0.045
N = 20
ax = Axis(fig[1:4, 10:13]; aspect=DataAspect(), spinewidth=4)
plot_z_quadrants!(dx*(N-1)/2/R/J, N, Rancic_backward, ax; box=true)

N = 24
ax = Axis(fig[6:9, 10:13]; aspect=DataAspect(), spinewidth=4)
plot_z_quadrants!(dx*(N-1)/2/R/J, N, MITgcm_backward, ax; box=true)
save(joinpath(@__DIR__, "grid_comparison.pdf"), fig)

# dx = 0.045
# N = 30
# L_km = dx * (N-1) / 2
# L = L_km / R / J
# fig, ax = plot_quadrants(L, N)
# save("mitgcm_map_zoom.pdf", fig)
