function plot1D(
    fname, var;
    zu=nothing, tu=nothing, norm=true, colormap=nothing, colorrange=nothing, save=false,
    save_fname=nothing, new_window=false,
)
    fp = HDF5.h5open(fname, "r")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "fields/t")
    F = HDF5.read(fp, "fields/" * string(var))
    HDF5.close(fp)

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

    isnothing(colorrange) ? colorrange = (minimum(F), maximum(F)) : nothing
    vmin, vmax = colorrange
    isnothing(vmin) ? vmin = minimum(F) : nothing
    isnothing(vmax) ? vmax = maximum(F) : nothing
    if vmin * vmax < 0   # diverging colormap
        isnothing(colormap) ? colormap = CMAPDIV : nothing
    else
        isnothing(colormap) ? colormap = CMAP : nothing
    end

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis(fig[1,1]; xlabel="t ($stu)", ylabel="z ($szu)")

    hm = mak.heatmap!(ax, t, z, transpose(F); colormap, colorrange)
    mak.Colorbar(fig[2,1], hm; vertical=false, label=string(var), flipaxis=false)

    if save
        if isnothing(save_fname)
            ext = splitext(fname)[end]
            save_fname = replace(fname, ext => ".png")
        end
        mak.save(save_fname, fig)
    end

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end
    return nothing
end


function plot1D_line(
    fname, var;
    zu=nothing, norm=true, vmin=nothing, vmax=nothing, save=false, save_fname=nothing,
    new_window=false,
)
    fp = HDF5.h5open(fname, "r")
    z = HDF5.read(fp, "z")
    F = HDF5.read(fp, string(var))
    HDF5.close(fp)

    zu = isnothing(zu) ? units(z) : zu
    szu = units_name_space(zu)

    @. z = z / zu

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    isnothing(vmin) ? vmin = minimum(F) : nothing
    isnothing(vmax) ? vmax = maximum(F) : nothing

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis(fig[1,1]; xlabel="z ($szu)", ylabel=string(var))

    mak.lines!(ax, z, F)

    if save
        if isnothing(save_fname)
            ext = splitext(fname)[end]
            save_fname = replace(fname, ext => ".png")
        end
        mak.save(save_fname, fig)
    end

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end
    return nothing
end
