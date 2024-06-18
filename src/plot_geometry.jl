function plot_geometry(fname; kwargs...)
    n = grid_dims(fname)
    if n == 1
        plot_geometry_1D(fname; kwargs...)
    elseif n == 2
        plot_geometry_2D(fname; kwargs...)
    elseif n == 3
        plot_geometry_3D(fname; kwargs...)
    end
    return nothing
end


# ******************************************************************************************
# 1D
# ******************************************************************************************
function plot_geometry_1D(fname; kwargs...)
    fp = HDF5.h5open(fname, "r")
    z = HDF5.read(fp, "z")
    geometry = HDF5.read(fp, "geometry")
    HDF5.close(fp)
    plot_geometry(z, geometry; kwargs...)
    return nothing
end


function plot_geometry(
    zin, geometry;
    zu=nothing, zlims=nothing, color=:gray15, new_window=false, fill=false, save=false,
    save_fname="out.png",
)
    iz1, iz2 = indices_of_limits(zin, zlims)
    z = zin[iz1:iz2]

    if typeof(geometry) <: Function
        F = [geometry(zi) ? 1 : 0 for zi=z]
    else
        F = [Int(geometry[iz]) for iz=iz1:iz2]
    end

    isnothing(zu) ? zu = units(z) : nothing
    szu = units_name_space(zu)

    z = z / zu

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis(fig[1,1]; xlabel="z ($szu)")
    mak.lines!(ax, z, F; color)
    if fill
        mak.fill_between!(ax, z, F, zero(F); color)
    end

    if save
        mak.save(save_fname, fig)
    end

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end
    return nothing
end


# ******************************************************************************************
# 2D
# ******************************************************************************************
function plot_geometry_2D(fname; kwargs...)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    z = HDF5.read(fp, "z")
    geometry = HDF5.read(fp, "geometry")
    HDF5.close(fp)
    plot_geometry(x, z, geometry; kwargs...)
    return nothing
end


function plot_geometry(
    xin, zin, geometry;
    xu=nothing, zu=nothing, xlims=nothing, zlims=nothing, aspect=1, colormap=nothing,
    colorrange=nothing, colorbar=true, new_window=false, save=false, save_fname="out.png",
)
    ix1, ix2 = indices_of_limits(xin, xlims)
    iz1, iz2 = indices_of_limits(zin, zlims)
    x = xin[ix1:ix2]
    z = zin[iz1:iz2]

    if typeof(geometry) <: Function
        F = [geometry(xi,zi) ? 1 : 0 for xi=x, zi=z]
    else
        F = [Int(geometry[ix,iz]) for ix=ix1:ix2, iz=iz1:iz2]
    end

    xu = isnothing(xu) ? units(x) : xu
    zu = isnothing(zu) ? units(z) : zu
    sxu = units_name_space(xu)
    szu = units_name_space(zu)

    x = x / xu
    z = z / zu

    nmat = maximum(F)   # number of materials

    colormap = isnothing(colormap) ? CMAP : colormap
    colormap = mak.to_colormap(colormap)
    colormap[1] = mak.RGBAf(0,0,0,0)

    colorrange = isnothing(colorrange) ? (0,nmat) : colorrange

    lmin = max(0, minimum(colorrange))
    lmax = min(nmat, maximum(colorrange)) + 1
    levels = lmin:lmax

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis(
        fig[1,1]; xlabel="x ($sxu)", ylabel="z ($szu)", aspect, xautolimitmargin=(0,0),
        yautolimitmargin=(0,0),
    )
    hm = mak.contourf!(ax, x, z, F; colormap, colorrange, levels)
    if colorbar
        mak.Colorbar(
            fig[2,1], hm;
            vertical=false, flipaxis=false, width=mak.Relative(3/4), ticks=lmin:lmax,
        )
    end

    if save
        mak.save(save_fname, fig)
    end

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end
    return nothing
end



# ******************************************************************************************
# 3D
# ******************************************************************************************
function plot_geometry_3D(fname; kwargs...)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    geometry = HDF5.read(fp, "geometry")
    HDF5.close(fp)
    plot_geometry(x, y, z, geometry; kwargs...)
    return nothing
end


function plot_geometry(
    xin, yin, zin, geometry;
    xu=nothing, yu=nothing, zu=nothing, xlims=nothing, ylims=nothing, zlims=nothing,
    aspect=(1,1,1), colormap=nothing, colorrange=nothing, colorbar=true, new_window=false,
    save=false, save_fname="out.png", lscene=true,
)
    ix1, ix2 = indices_of_limits(xin, xlims)
    iy1, iy2 = indices_of_limits(yin, ylims)
    iz1, iz2 = indices_of_limits(zin, zlims)
    x = xin[ix1:ix2]
    y = yin[iy1:iy2]
    z = zin[iz1:iz2]

    if typeof(geometry) <: Function
        F = [geometry(xi,yi,zi) ? 1 : 0 for xi=x, yi=y, zi=z]
    else
        F = [Int(geometry[ix,iy,iz]) for ix=ix1:ix2, iy=iy1:iy2, iz=iz1:iz2]
    end

    xu = isnothing(xu) ? units(x) : xu
    yu = isnothing(yu) ? units(y) : yu
    zu = isnothing(zu) ? units(z) : zu
    sxu = units_name_space(xu)
    syu = units_name_space(yu)
    szu = units_name_space(zu)

    x = x / xu
    y = y / yu
    z = z / zu

    nmat = maximum(F)   # number of materials

    colormap = isnothing(colormap) ? CMAP : colormap
    colorrange = isnothing(colorrange) ? (0,nmat) : colorrange

    lmin = max(0, minimum(colorrange)) + 1
    lmax = min(nmat, maximum(colorrange)) + 1
    levels = lmin:lmax

    fig = mak.Figure(size=(950,992))
    if lscene
        scale = scene_scale(x, y, z, aspect)

        ax = mak.LScene(fig[1,1])
        mak.scale!(ax.scene, scale)
        axis = ax.scene[mak.OldAxis]
        axis.names[:axisnames] = ("x ($sxu)", "y ($syu)", "z ($szu)")
        axis.names[:fontsize] = 3
        axis.ticks[:fontsize] = 3
    else
        ax = mak.Axis3(
            fig[1,1]; xlabel="x ($sxu)", ylabel="y ($syu)", zlabel="z ($szu)", aspect,
            xautolimitmargin=(0,0), yautolimitmargin=(0,0), zautolimitmargin=(0,0),
            perspectiveness=0,
        )
    end

    xrange = (x[1], x[end])
    yrange = (y[1], y[end])
    zrange = (z[1], z[end])

    img = mak.contour!(ax, xrange, yrange, zrange, F; colormap, colorrange, levels)
    if colorbar
        mak.Colorbar(fig[1,2], img; height=mak.Relative(3/4), ticks=lmin:lmax)
    end

    if save
        mak.save(save_fname, fig)
    end

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end
    return nothing
end
