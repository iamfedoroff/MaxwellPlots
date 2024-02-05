function plot_viewpoints(
    fname, var, viewpoints; new_window=false, norm=true, xlims=nothing, ylims=nothing,
    zu=nothing, tu=nothing, Fu=1, save=false, save_fname="image.png",
)
    Np = length(viewpoints)

    fp = HDF5.h5open(fname, "r")

    if ! ("viewpoints" in keys(fp))
        error("No viewpoints in $fname")
    end

    group = fp["viewpoints"]

    for n=1:Np
        if ! (string(viewpoints[n]) in keys(group))
            error("Viewpoint $(viewpoints[n]) not in $fname")
        end
    end

    z = HDF5.read(fp, "z")
    t = HDF5.read(group, "t")

    isnothing(zu) ? zu = units(z) : nothing
    isnothing(tu) ? tu = units(t) : nothing
    szu = units_name_space(zu)
    stu = units_name_time(tu)

    Enorm = 1

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis(
        fig[1,1]; xlabel="t ($stu)", ylabel=string(var), subtitle=basename(fname),
        limits=(xlims,ylims),
    )

    for (n, pt) in enumerate(viewpoints)
        point = HDF5.read(group, "$pt/point")
        Et = HDF5.read(group, "$pt/"*string(var))

        if norm
            if n == 1
                Enorm = maximum(Et)
            end
            @. Et = Et / Enorm
        else
            @. Et = Et / Fu
        end

        coords = Tuple([round(coord/zu; digits=2) for coord in point])
        label = "$coords $szu"

        mak.lines!(ax, t/tu, Et; label)
    end

    mak.axislegend(ax; framevisible=false)

    HDF5.close(fp)

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end

    if save
        mak.save(save_fname, fig)
    end

    return nothing
end


function plot_viewpoints_spectrum(
    fname, var, viewpoints; new_window=false, norm=true, xlims=nothing, ylims=nothing,
    zu=nothing, wu=1, Fu=1, save=false, yscale=log10,
)
    Np = length(viewpoints)

    fp = HDF5.h5open(fname, "r")

    if ! ("viewpoints" in keys(fp))
        error("No viewpoints in $fname")
    end

    group = fp["viewpoints"]

    for n=1:Np
        if ! (string(viewpoints[n]) in keys(group))
            error("Viewpoint $(viewpoints[n]) not in $fname")
        end
    end

    z = HDF5.read(fp, "z")
    t = HDF5.read(group, "t")

    Nt = length(t)
    dt = t[2] - t[1]
    w = 2*pi * FFTW.rfftfreq(Nt, 1/dt)

    isnothing(zu) ? zu = units(z) : nothing
    szu = units_name_space(zu)

    Snorm = 1

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis(
        fig[1,1];
        xlabel="Frequency w", ylabel="Spectrum of " * string(var), subtitle=basename(fname),
        limits=(xlims,ylims), yscale,
    )

    for (n, pt) in enumerate(viewpoints)
        point = HDF5.read(group, "$pt/point")
        Et = HDF5.read(group, "$pt/"*string(var))

        Ew = FFTW.rfft(Et)
        Sw = @. abs2(Ew)
        if norm
            if n == 1
                Snorm = maximum(Sw)
            end
            @. Sw = Sw / Snorm
        else
            @. Sw = Sw / Fu
        end

        if yscale == log10
            Sw = replace(Sw, 0=>eps(eltype(Sw)))   # to avoid log10(0)
        end

        coords = Tuple([round(coord/zu; digits=2) for coord in point])
        label = "$coords $szu"

        mak.lines!(ax, w/wu, Sw; label)
    end

    mak.axislegend(ax; framevisible=false)

    HDF5.close(fp)

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end

    if save
        ext = splitext(fname)[end]
        fname_fig = replace(fname, ext => ".png")
        mak.save(fname_fig, fig)
    end

    return nothing
end


function plot_viewpoints_polarization(
    fname, var1, var2, viewpoints;
    zu=nothing, tu=nothing, tlims=nothing, Flims=(-1,1), save=false, save_fname="image.png",
    new_window=false,
)
    Np = length(viewpoints)

    fp = HDF5.h5open(fname, "r")

    if ! ("viewpoints" in keys(fp))
        error("No viewpoints in $fname")
    end

    group = fp["viewpoints"]

    for n=1:Np
        if ! (string(viewpoints[n]) in keys(group))
            error("Viewpoint $(viewpoints[n]) not in $fname")
        end
    end

    z = HDF5.read(fp, "z")
    t = HDF5.read(group, "t")

    isnothing(zu) ? zu = units(z) : nothing
    isnothing(tu) ? tu = units(t) : nothing
    szu = units_name_space(zu)
    stu = units_name_time(tu)

    @. t = t / tu

    if isnothing(tlims)
        tmin, tmax = extrema(t)
    else
        tmin, tmax = tlims
    end
    mask = @. (t >= tmin) && (t <= tmax)
    t = t[mask]

    # read data to find the normalization factor:
    Enorm = 0
    for pt in viewpoints
        Ex = HDF5.read(group, string(pt) * "/" * string(var1))
        Ey = HDF5.read(group, string(pt) * "/" * string(var2))
        Exmin, Exmax = extrema(Ex)
        Eymin, Eymax = extrema(Ey)
        Emin = sqrt(Exmin^2 + Eymin^2)
        Emax = sqrt(Exmax^2 + Eymax^2)
        Enorm = max(Enorm, Emin, Emax)
    end
    @show Enorm

    xmin, xmax = tmin, tmax
    ymin, ymax = Flims
    zmin, zmax = Flims

    tones = ones(length(t))

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis3(
        fig[1,1]; perspectiveness=0,
        xlabel="t ($stu)", ylabel=string(var2), zlabel=string(var1),
    )
    mak.xlims!(ax, (xmin,xmax))
    mak.ylims!(ax, (ymin,ymax))
    mak.zlims!(ax, (zmin,zmax))
    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end

    # read data, normalize, and plot:
    for (n, pt) in enumerate(viewpoints)
        point = HDF5.read(group, "$pt/point")
        Ex = HDF5.read(group, string(pt) * "/" * string(var1))
        Ey = HDF5.read(group, string(pt) * "/" * string(var2))

        @. Ex = Ex / Enorm
        @. Ey = Ey / Enorm

        Ex = Ex[mask]
        Ey = Ey[mask]

        coords = Tuple([round(coord/zu; digits=2) for coord in point])
        label = "$coords $szu"

        mak.lines!(ax, t, Ey, Ex; label)
        mak.lines!(ax, t, zmax * tones, Ex; color=mak.Cycled(n), linewidth=0.5)
        mak.lines!(ax, t, Ey, ymin * tones; color=mak.Cycled(n), linewidth=0.5)
        mak.lines!(ax, xmax * tones, Ey, Ex; color=mak.Cycled(n), linewidth=0.5)
    end

    mak.axislegend(ax; framevisible=false)

    HDF5.close(fp)

    if save
        mak.save(save_fname, fig)
    end

    return nothing
end
