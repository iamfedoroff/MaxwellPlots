function plot_geometry(fname; kw_args...)
    fp = HDF5.h5open(fname, "r")
    if "y" in keys(fp)
        # 3D:
        dims = 3
        x = HDF5.read(fp, "x")
        y = HDF5.read(fp, "y")
        z = HDF5.read(fp, "z")
    elseif "x" in keys(fp)
        # 2D:
        dims = 2
        x = HDF5.read(fp, "x")
        z = HDF5.read(fp, "z")
    else
        # 1D:
        dims = 1
        z = HDF5.read(fp, "z")
    end
    geometry = HDF5.read(fp, "geometry")
    HDF5.close(fp)

    if dims == 3
        # 3D:
        plot_geometry(x, y, z, geometry; kw_args...)
    elseif dims == 2
        # 2D:
        plot_geometry(x, z, geometry; kw_args...)
    else
        # 1D:
        plot_geometry(z, geometry; kw_args...)
    end
    return nothing
end


function plot_geometry(
    z, geometry; zu=nothing, new_window=false, color=:gray15, fill=false, save=false,
    save_fname="out.png",
)
    if typeof(geometry) <: Function
        F = [geometry(zi) ? 1 : 0 for zi=z]
    else
        F = geometry
    end

    isnothing(zu) ? zu = units(z) : nothing
    szu = units_name_space(zu)

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis(fig[1,1]; xlabel="z ($szu)", ylabel="Material geometry")
    mak.lines!(ax, z/zu, F; color)
    if fill
        mak.fill_between!(ax, z/zu, F, zero(F); color)
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


function plot_geometry(
    x, z, geometry; xu=nothing, zu=nothing, cmap=CMAP, aspect=1, new_window=false,
    save=false, save_fname="out.png",
)
    if typeof(geometry) <: Function
        F = [geometry(xi,zi) ? 1 : 0 for xi=x, zi=z]
    else
        F = geometry
    end

    xu = isnothing(xu) ? units(x) : xu
    zu = isnothing(zu) ? units(z) : zu
    sxu = units_name_space(xu)
    szu = units_name_space(zu)

    xx = x / xu
    zz = z / zu

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis(
        fig[1,1]; xlabel="x ($sxu)", ylabel="z ($szu)", title="Material geometry", aspect,
    )
    hm = mak.heatmap!(
        ax, xx, zz, F; colormap=cmap, colorrange=(0,1),
    )
    # mak.Colorbar(fig[2,1], hm; vertical=false, label="geometry")

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


function plot_geometry(
    x, y, z, geometry;
    xu=nothing, yu=nothing, zu=nothing, xlims=nothing, ylims=nothing, zlims=nothing,
    colormap=CMAP, aspect=(1,1,1), new_window=false, algorithm=:iso, isovalue=1,
    absorption=1, save=false, save_fname="out.png",
)
    if typeof(geometry) <: Function
        F = [geometry(xi,yi,zi) ? 1 : 0 for xi=x, yi=y, zi=z]
    else
        F = geometry
    end

    x, y, z, F = apply_limits(x, y, z, F; xlims, ylims, zlims)

    xu = isnothing(xu) ? units(x) : xu
    yu = isnothing(yu) ? units(y) : yu
    zu = isnothing(zu) ? units(z) : zu
    sxu = units_name_space(xu)
    syu = units_name_space(yu)
    szu = units_name_space(zu)

    xx = x / xu
    yy = y / yu
    zz = z / zu

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis3(
        fig[1,1];
        xlabel="x ($sxu)", ylabel="y ($syu)", zlabel="z ($szu)",
        title="Material geometry", aspect, perspectiveness=0,
    )
    mak.limits!(ax, xx[1], xx[end], yy[1], yy[end], zz[1], zz[end])
    img = mak.volume!(ax, xx, yy, zz, F; colormap, algorithm, isovalue, absorption)
    # mak.Colorbar(fig[1,2], img; label="geometry")

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
