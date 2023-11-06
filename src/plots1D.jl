function plot1D(
    fname, var;
    zu=1, tu=1, norm=true, colormap=nothing, colorrange=nothing, save=false,
    save_fname=nothing, new_window=false,
)
    fp = HDF5.h5open(fname, "r")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "fields/t")
    F = HDF5.read(fp, "fields/" * string(var))
    HDF5.close(fp)

    @. z = z / zu
    @. t = t / tu
    szu = space_units_string(zu)
    stu = time_units_string(tu)

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

    fig = mak.Figure(resolution=(950,992), fontsize=14)
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


function inspect1D(fname, var; zu=1, tu=1, vmin=-1, vmax=1, norm=true, new_window=false)
    fp = HDF5.h5open(fname, "r")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "fields/t")
    F = HDF5.read(fp, "fields/" * string(var))
    HDF5.close(fp)

    @. z = z / zu
    @. t = t / tu
    szu = space_units_string(zu)
    stu = time_units_string(tu)

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    fig = mak.Figure(resolution=(950,992), fontsize=14)
    ax = mak.Axis(fig[1,1]; xlabel="z ($szu)", ylabel=string(var))

    it = 1
    line = mak.lines!(ax, z, F[:,it])
    mak.ylims!(ax, (vmin,vmax))
    ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)

    sg = mak.SliderGrid(fig[2,1], (label="Time", range=1:length(t), startvalue=1))
    mak.on(sg.sliders[1].value) do it
        line[2] = F[:,it]
        ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
    end

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end
    return nothing
end
