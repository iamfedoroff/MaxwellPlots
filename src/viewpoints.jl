function plot_viewpoints(
    fname, svar, viewpoints; new_window=false, norm=false, xlims=nothing, ylims=nothing,
    zu=1, tu=1, Fu=1, save=false,
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

    t = HDF5.read(group, "t")

    szu = space_units_string(zu)
    stu = time_units_string(tu)

    Enorm = 1

    fig = mak.Figure()
    ax = mak.Axis(
        fig[1,1]; xlabel="t ($stu)", ylabel=svar, subtitle=basename(fname),
        limits=(xlims,ylims),
    )

    for (n, pt) in enumerate(viewpoints)
        point = HDF5.read(group, "$pt/point")
        Et = HDF5.read(group, "$pt/"*svar)

        if norm
            # if n == 1
            #     Enorm = maximum(Et)
            # end
            Enorm = maximum(Et)
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
        ext = splitext(fname)[end]
        fname_fig = replace(fname, ext => ".png")
        mak.save(fname_fig, fig)
    end

    return nothing
end


function plot_viewpoints_spectrum(
    fname, svar, viewpoints; new_window=false, norm=true, xlims=nothing, ylims=nothing,
    zu=1, wu=1, Fu=1, save=false, yscale=log10,
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

    t = HDF5.read(group, "t")
    Nt = length(t)
    dt = t[2] - t[1]
    w = 2*pi * FFTW.rfftfreq(Nt, 1/dt)

    szu = space_units_string(zu)

    Snorm = 1

    fig = mak.Figure()
    ax = mak.Axis(
        fig[1,1];
        xlabel="Frequency w", ylabel="Spectrum of " * svar, subtitle=basename(fname),
        limits=(xlims,ylims), yscale,
    )

    for (n, pt) in enumerate(viewpoints)
        point = HDF5.read(group, "$pt/point")
        Et = HDF5.read(group, "$pt/"*svar)

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
