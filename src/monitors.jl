function inds2array(inds)
    return [CartesianIndex(ind[1]...) for ind in inds]
end


function plot_monitors(
    fname, var, monitors; new_window=false, norm=true, xlims=nothing, ylims=nothing,
    zu=nothing, tu=nothing, Fu=1, save=false, save_fname="out.png",
)
    Nm = length(monitors)

    fp = HDF5.h5open(fname, "r")

    if ! ("monitors" in keys(fp))
        error("No monitors in $fname")
    end

    group = fp["monitors"]

    for n=1:Nm
        if ! (string(monitors[n]) in keys(group))
            error("Monitor $(monitors[n]) is not in $fname")
        end
    end

    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "t")

    zu = isnothing(zu) ? units(z) : zu
    tu = isnothing(tu) ? units(t) : tu
    szu = units_name_space(zu)
    stu = units_name_time(tu)

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis(
        fig[1,1]; xlabel="t ($stu)", ylabel=string(var), subtitle=basename(fname),
        limits=(xlims,ylims),
    )

    Fnorm = 1

    for (n, mon) in enumerate(monitors)
        issum = HDF5.read(group, "$mon/issum")
        # inds = HDF5.read(group, "$mon/inds")
        # ci = inds2array(inds)

        Ft = HDF5.read(group, "$mon/"*string(var))
        Ft = issum ? Ft : Ft[1,:]

        @show extrema(Ft)
        if norm
            if n == 1
                Fnorm = maximum(Ft)
            end
            @. Ft = Ft / Fnorm
        else
            @. Ft = Ft / Fu
        end

        mak.lines!(ax, t/tu, Ft; label="$n")
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


function plot_monitors_spectrum(
    fname, var, monitors; new_window=false, norm=true, xlims=nothing, ylims=nothing,
    zu=nothing, wu=1, Fu=1, save=false, save_fname="out.png", yscale=log10,
)
    Nm = length(monitors)

    fp = HDF5.h5open(fname, "r")

    if ! ("monitors" in keys(fp))
        error("No monitors in $fname")
    end

    group = fp["monitors"]

    for n=1:Nm
        if ! (string(monitors[n]) in keys(group))
            error("Monitor $(monitors[n]) is not in $fname")
        end
    end

    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "t")

    Nt = length(t)
    dt = t[2] - t[1]
    w = 2*pi * FFTW.rfftfreq(Nt, 1/dt)

    zu = isnothing(zu) ? units(z) : zu
    szu = units_name_space(zu)

    Snorm = 1

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis(
        fig[1,1];
        xlabel="Frequency w", ylabel="Spectrum of " * string(var), subtitle=basename(fname),
        limits=(xlims,ylims), yscale,
    )

    for (n, mon) in enumerate(monitors)
        issum = HDF5.read(group, "$mon/issum")
        # inds = HDF5.read(group, "$mon/inds")
        # ci = inds2array(inds)

        Et = HDF5.read(group, "$mon/"*string(var))
        Et = issum ? Et : Et[1,:]

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

        mak.lines!(ax, w/wu, Sw; label="$n")
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


function plot_monitors_polarization(
    fname, var1, var2, monitors;
    zu=nothing, tu=nothing, tlims=nothing, Flims=(-1,1), save=false, save_fname="out.png",
    new_window=false,
)
    Nm = length(monitors)

    fp = HDF5.h5open(fname, "r")

    if ! ("monitors" in keys(fp))
        error("No monitors in $fname")
    end

    group = fp["monitors"]

    for n=1:Nm
        if ! (string(monitors[n]) in keys(group))
            error("Monitor $(monitors[n]) is not in $fname")
        end
    end

    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "t")

    zu = isnothing(zu) ? units(z) : zu
    tu = isnothing(tu) ? units(t) : tu
    szu = units_name_space(zu)
    stu = units_name_time(tu)

    if isnothing(tlims)
        tmin, tmax = extrema(t)
    else
        tmin, tmax = tlims
    end
    mask = @. (t >= tmin) && (t <= tmax)
    t = t[mask]

    # read data to find the normalization factor:
    Enorm = 0
    for mon in monitors
        Ex = HDF5.read(group, string(mon) * "/" * string(var1))
        Ey = HDF5.read(group, string(mon) * "/" * string(var2))
        Exmin, Exmax = extrema(Ex)
        Eymin, Eymax = extrema(Ey)
        Emin = sqrt(Exmin^2 + Eymin^2)
        Emax = sqrt(Exmax^2 + Eymax^2)
        Enorm = max(Enorm, Emin, Emax)
    end
    @show Enorm

    xmin, xmax = tmin/tu, tmax/tu
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

    # read data, normalize, and plot:
    for (n, mon) in enumerate(monitors)
        issum = HDF5.read(group, "$mon/issum")
        # inds = HDF5.read(group, "$mon/inds")
        # ci = inds2array(inds)
        Ex = HDF5.read(group, string(mon) * "/" * string(var1))
        Ey = HDF5.read(group, string(mon) * "/" * string(var2))

        Ex = issum ? Ex : Ex[1,:]
        Ey = issum ? Ey : Ey[1,:]
        @. Ex = Ex / Enorm
        @. Ey = Ey / Enorm
        Ex = Ex[mask]
        Ey = Ey[mask]

        mak.lines!(ax, t/tu, Ey, Ex; label="$n")
        mak.lines!(ax, t/tu, zmax * tones, Ex; color=mak.Cycled(n), linewidth=0.5)
        mak.lines!(ax, t/tu, Ey, ymin * tones; color=mak.Cycled(n), linewidth=0.5)
        mak.lines!(ax, xmax * tones, Ey, Ex; color=mak.Cycled(n), linewidth=0.5)
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
