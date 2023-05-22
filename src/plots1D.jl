function plot1D(
    fname, svar;
    zu=1, tu=1, vmin=-1, vmax=1, norm=false, cmap=:seismic, new_window=false, save=false,
)
    fp = HDF5.h5open(fname, "r")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "t")
    F = HDF5.read(fp, svar)
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
    if new_window
        mak.display(mak.Screen(), fig)
    end
    ax = mak.Axis(fig[1,1]; xlabel="t ($stu)", ylabel="z ($szu)")
    mak.display(fig)

    hm = mak.heatmap!(
        ax, t, z, transpose(F); colormap=cmap, colorrange=(vmin,vmax),
    )
    mak.Colorbar(fig[2,1], hm; vertical=false, label=svar, flipaxis=false)

    if save
        ext = splitext(fname)[end]
        fname_fig = replace(fname, ext => ".png")
        mak.save(fname_fig, fig)
    end
    return nothing
end


function inspect1D(fname, svar; zu=1, tu=1, vmin=-1, vmax=1, norm=false)
    fp = HDF5.h5open(fname, "r")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "t")
    F = HDF5.read(fp, svar)
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
    ax = mak.Axis(fig[1,1]; xlabel="z ($szu)", ylabel=svar)
    mak.display(fig)

    it = 1
    line = mak.lines!(ax, z, F[:,it])
    mak.ylims!(ax, (vmin,vmax))
    ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)

    sg = mak.SliderGrid(fig[2,1], (label="Time", range=1:length(t), startvalue=1))
    mak.on(sg.sliders[1].value) do it
        line[2] = F[:,it]
        ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
    end
    return nothing
end
