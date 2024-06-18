function inspect(fname, var; kwargs...)
    n = grid_dims(fname)
    if n == 1
        inspect1D(fname, var; kwargs...)
    elseif n == 2
        inspect2D(fname, var; kwargs...)
    elseif n == 3
        inspect3D(fname, var; kwargs...)
    end
    return nothing
end


# ******************************************************************************************
# 1D
# ******************************************************************************************
function inspect1D(fname, var; kwargs...)
    fp = HDF5.h5open(fname, "r")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "fields/t")
    F = HDF5.read(fp, "fields/" * string(var))
    HDF5.close(fp)
    inspect(z, t, F; label=string(var), kwargs...)
    return nothing
end


function inspect(
    zin, tin, Fin;
    zu=nothing, tu=nothing, zlims=nothing, tlims=nothing, norm=true, new_window=false,
    movie=false, movie_fname="out.mp4", movie_framerate=10, label="",
)
    z, t, F = apply_limits(zin, tin, Fin; xlims=zlims, ylims=tlims)

    zu = isnothing(zu) ? units(z) : zu
    tu = isnothing(tu) ? units(t) : tu
    szu = units_name_space(zu)
    stu = units_name_time(tu)

    @. z = z / zu
    @. t = t / tu

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis(fig[1,1]; xlabel="z ($szu)", ylabel=label)

    it = 1
    line = mak.lines!(ax, z, F[:,it])
    ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)

    if movie
        mak.record(fig, movie_fname, 1:length(t); framerate=movie_framerate) do it
            line[2] = F[:,it]
            ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end
    else
        sg = mak.SliderGrid(fig[2,1], (label="Time", range=1:length(t), startvalue=1))
        mak.on(sg.sliders[1].value) do it
            line[2] = F[:,it]
            ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end

        slider_keyboard_control(fig, sg.sliders[1])

        if new_window
            mak.display(mak.Screen(), fig)
        else
            mak.display(fig)
        end
    end

    return nothing
end


# ******************************************************************************************
# 2D
# ******************************************************************************************
function inspect2D(fname, var; kwargs...)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "fields/t")
    if var in (:Ex, :Ez, :Hy, :rho)
        F = HDF5.read(fp, "fields/" * string(var))
    elseif var == :poynting
        Hy = HDF5.read(fp, "fields/Hy")
        Ex = HDF5.read(fp, "fields/Ex")
        Ez = HDF5.read(fp, "fields/Ez")
        F = poynting(Hy, Ex, Ez)
    else
        error("Wrong input varible " * string(var))
    end
    HDF5.close(fp)

    inspect(x, z, t, F; label=string(var), kwargs...)
    return nothing
end


function inspect(
    xin, zin, tin, Fin;
    xu=nothing, zu=nothing, tu=nothing, xlims=nothing, zlims=nothing, tlims=nothing,
    norm=true, colormap=nothing, colorrange=nothing, colorbar=true, aspect=1,
    new_window=false, movie=false, movie_fname="out.mp4", movie_framerate=10, label="",
)
    x, z, t, F = apply_limits(xin, zin, tin, Fin; xlims=xlims, ylims=zlims, zlims=tlims)

    xu = isnothing(xu) ? units(x) : xu
    zu = isnothing(zu) ? units(z) : zu
    tu = isnothing(tu) ? units(t) : tu
    sxu = units_name_space(xu)
    szu = units_name_space(zu)
    stu = units_name_time(tu)

    @. x = x / xu
    @. z = z / zu
    @. t = t / tu

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    colormap = isnothing(colormap) ? default_colormap(F) : colormap
    colorrange = isnothing(colorrange) ? default_colorrange(F) : colorrange

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="z ($szu)", aspect)

    it = 1
    hm = mak.heatmap!(ax, x, z, F[:,:,it]; colormap, colorrange)
    if colorbar
        mak.Colorbar(
            fig[2,1], hm; label, vertical=false, flipaxis=false, width=mak.Relative(3/4),
        )
    end
    ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)

    if movie
        mak.record(fig, movie_fname, 1:length(t); framerate=movie_framerate) do it
            hm[3] = F[:,:,it]
            ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end
    else
        sg = mak.SliderGrid(fig[3,1], (label="Time", range=1:length(t), startvalue=1))
        mak.on(sg.sliders[1].value) do it
            hm[3] = F[:,:,it]
            ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end

        slider_keyboard_control(fig, sg.sliders[1])

        if new_window
            mak.display(mak.Screen(), fig)
        else
            mak.display(fig)
        end
    end

    return nothing
end


# ******************************************************************************************
# 3D
# ******************************************************************************************
function inspect3D(fname, var; kwargs...)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "fields/t")
    if string(var) in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "rho")
        F = HDF5.read(fp, "fields/" * string(var))
    elseif string(var) == "poynting"
        Hx = HDF5.read(fp, "fields/Hx")
        Hy = HDF5.read(fp, "fields/Hy")
        Hz = HDF5.read(fp, "fields/Hz")
        Ex = HDF5.read(fp, "fields/Ex")
        Ey = HDF5.read(fp, "fields/Ey")
        Ez = HDF5.read(fp, "fields/Ez")
        F = poynting(Hx, Hy, Hz, Ex, Ey, Ez)
    elseif string(var) == "divE"
        Ex = HDF5.read(fp, "fields/Ex")
        Ey = HDF5.read(fp, "fields/Ey")
        Ez = HDF5.read(fp, "fields/Ez")
        F = divergence(x, y, z, Ex, Ey, Ez)
    elseif string(var) == "divH"
        Hx = HDF5.read(fp, "fields/Hx")
        Hy = HDF5.read(fp, "fields/Hy")
        Hz = HDF5.read(fp, "fields/Hz")
        F = divergence(x, y, z, Hx, Hy, Hz)
    else
        error("Wrong input varible " * string(var))
    end
    HDF5.close(fp)

    inspect(x, y, z, t, F; label=string(var), kwargs...)
    return nothing
end


function inspect(
    xin, yin, zin, tin, Fin;
    xu=nothing, yu=nothing, zu=nothing, tu=nothing, xlims=nothing, ylims=nothing,
    zlims=nothing, tlims=nothing, norm=true, colormap=nothing, colorrange=nothing,
    colorbar=true, aspect=(1,1,1), new_window=false, movie=false, movie_fname="out.mp4",
    movie_framerate=10, lscene=true, label="",
)
    x, y, z, t, F = apply_limits(xin, yin, zin, tin, Fin; xlims, ylims, zlims, tlims)

    xu = isnothing(xu) ? units(x) : xu
    yu = isnothing(yu) ? units(y) : yu
    zu = isnothing(zu) ? units(z) : zu
    tu = isnothing(tu) ? units(t) : tu
    sxu = units_name_space(xu)
    syu = units_name_space(yu)
    szu = units_name_space(zu)
    stu = units_name_time(tu)

    @. x = x / xu
    @. y = y / yu
    @. z = z / zu
    @. t = t / tu

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    colormap = isnothing(colormap) ? default_colormap(F) : colormap
    colorrange = isnothing(colorrange) ? default_colorrange(F) : colorrange

    colormap = mak.to_colormap(colormap)
    if prod(colorrange) < 0
        Nmid = halfint(length(colormap))
        @. colormap[Nmid-1:Nmid+1] = mak.RGBAf(0,0,0,0)
    else
        colormap[1] = mak.RGBAf(0,0,0,0)
    end

    fig = mak.Figure(size=(950,992))
    if lscene
        scale = scene_scale(x, y, z, aspect)

        ax = mak.LScene(fig[1,1])
        mak.scale!(ax.scene, scale)
        axis = ax.scene[mak.OldAxis]
        axis.names[:axisnames] = ("x ($sxu)", "y ($syu)", "z ($szu)")
        axis.names[:fontsize] = 3
        axis.ticks[:fontsize] = 3
        title = mak.Label(fig[1,1,mak.Top()], "").text

        # Does not work properly. Clips the data.
        # cam = mak.cameracontrols(ax)
        # cam.settings.projectiontype[] = Makie.Orthographic
        # mak.update_cam!(ax.scene, cam)
    else
        ax = mak.Axis3(
            fig[1,1]; aspect, perspectiveness=0, xlabel="x ($sxu)", ylabel="y ($syu)",
            zlabel="z ($szu)",
        )
        title = ax.title
    end

    xrange = (x[1], x[end])
    yrange = (y[1], y[end])
    zrange = (z[1], z[end])
    it = 1

    img = mak.volume!(
        ax, xrange, yrange, zrange, F[:,:,:,it]; colormap, colorrange,
        algorithm=:absorption, absorption=4f0,
        # algorithm=:iso, isovalue=0.5*vmax,
    )
    # img = mak.contour!(
    #     ax, xrange, yrange, zrange, F[:,:,:,it];
    #     levels=[0.1*vmin,0.1*vmax], colormap=cmap, colorrange=(vmin,vmax),
    #     alpha=1,
    # )
    if colorbar
        mak.Colorbar(fig[1,2], img; label, height=mak.Relative(0.5))
    end
    title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)

    if movie
        mak.record(fig, movie_fname, 1:length(t); framerate=movie_framerate) do it
            img[4] = F[:,:,:,it]
            title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end
    else
        sg = mak.SliderGrid(fig[2,1], (label="Time", range=1:length(t), startvalue=1))
        mak.on(sg.sliders[1].value) do it
            img[4] = F[:,:,:,it]
            title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end

        slider_keyboard_control(fig, sg.sliders[1])

        if new_window
            mak.display(mak.Screen(), fig)
        else
            mak.display(fig)
        end
    end

    return nothing
end
