function plot2D(
    fname, var, t0;
    xu=1, zu=1, tu=1, norm=true, colormap=nothing, colorrange=nothing, aspect=1,
    new_window=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "fields/t")
    F = HDF5.read(fp, "fields/" * string(var))
    HDF5.close(fp)

    @. x = x / xu
    @. z = z / zu
    @. t = t / tu
    sxu = space_units_string(xu)
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

    it = argmin(abs.(t .- t0))

    fig = mak.Figure(resolution=(950,992), fontsize=14)
    ax = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="z ($szu)", aspect)

    hm = mak.heatmap!(
        ax, x, z, F[:,:,it]; colormap, colorrange=(vmin,vmax),
    )
    mak.Colorbar(fig[2,1], hm; vertical=false, label=string(var), flipaxis=false)
    ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end
    return nothing
end


function inspect2D(
    fname, var;
    xu=1, zu=1, tu=1, xlims=(nothing,nothing), zlims=(nothing,nothing), norm=true,
    colormap=nothing, colorrange=nothing, aspect=1, movie=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "fields/t")
    if var in (:Ex, :Ez, :Hy)
        F = HDF5.read(fp, "fields/" * string(var))
        isnothing(colormap) ? colormap = CMAPDIV : nothing
    elseif var == :poynting
        Hy = HDF5.read(fp, "fields/Hy")
        Ex = HDF5.read(fp, "fields/Ex")
        Ez = HDF5.read(fp, "fields/Ez")
        F = poynting(Hy, Ex, Ez)
        isnothing(colormap) ? colormap = CMAP : nothing
    else
        error("Wrong input varible " * string(var))
    end
    HDF5.close(fp)

    @. x = x / xu
    @. z = z / zu
    @. t = t / tu
    sxu = space_units_string(xu)
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

    isnothing(xlims[1]) ? xmin=x[1] : xmin=xlims[1]
    isnothing(xlims[2]) ? xmax=x[end] : xmax=xlims[2]
    isnothing(zlims[1]) ? zmin=z[1] : zmin=zlims[1]
    isnothing(zlims[2]) ? zmax=z[end] : zmax=zlims[2]

    fig = mak.Figure(resolution=(950,992), fontsize=14)
    ax = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="z ($szu)", aspect)
    mak.xlims!(ax, (xmin,xmax))
    mak.ylims!(ax, (zmin,zmax))
    mak.display(fig)

    it = 1
    hm = mak.heatmap!(ax, x, z, F[:,:,it]; colormap, colorrange)
    mak.Colorbar(fig[2,1], hm; vertical=false, label=string(var), flipaxis=false)
    ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)

    if movie
        ext = splitext(fname)[end]
        fname_movie = replace(fname, ext => ".mp4")
        mak.record(fig, fname_movie, 1:length(t); framerate=12) do it
            hm[3] = F[:,:,it]
            ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end
    else
        sg = mak.SliderGrid(fig[3,1], (label="Time", range=1:length(t), startvalue=1))
        mak.on(sg.sliders[1].value) do it
            hm[3] = F[:,:,it]
            ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end
    end
    return nothing
end


function inspect2D_xsec(
    fname, var, x0, z0; xu=1, zu=1, tu=1, vmin=-1, vmax=1, norm=true,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "fields/t")
    F = HDF5.read(fp, "fields/" * string(var))
    HDF5.close(fp)

    @. x = x / xu
    @. z = z / zu
    @. t = t / tu
    sxu = space_units_string(xu)
    szu = space_units_string(zu)
    stu = time_units_string(tu)

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    ix0 = argmin(abs.(x .- x0))
    iz0 = argmin(abs.(z .- z0))

    fig = mak.Figure(resolution=(950,992), fontsize=14)
    ax1 = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel=string(var))
    ax2 = mak.Axis(fig[2,1]; xlabel="z ($szu)", ylabel=string(var))
    mak.display(fig)

    it = 1
    line1 = mak.lines!(ax1, x, F[:,iz0,it])
    line2 = mak.lines!(ax2, z, F[ix0,:,it])
    mak.ylims!(ax1, (vmin,vmax))
    mak.ylims!(ax2, (vmin,vmax))
    ax1.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)

    sg = mak.SliderGrid(fig[3,1], (label="Time", range=1:length(t), startvalue=1))
    mak.on(sg.sliders[1].value) do it
        line1[2] = F[:,iz0,it]
        line2[2] = F[ix0,:,it]
        ax1.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
    end
    return nothing
end


function plot2D_poynting_averaged(
    fname; xu=1, zu=1, vmin=0, vmax=1, norm=false, norm_point=nothing, aspect=1,
    xlims=(nothing,nothing), zlims=(nothing,nothing),
    cmap=CMAP,
    new_window=false, save=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    z = HDF5.read(fp, "z")
    F = HDF5.read(fp, "Sa")
    HDF5.close(fp)

    @. x = x / xu
    @. z = z / zu
    sxu = space_units_string(xu)
    szu = space_units_string(zu)

    @show extrema(F)
    if norm
        if isnothing(norm_point)
            F .= F ./ maximum(F)
        else
            xn, zn = norm_point
            ixn = argmin(abs.(x .- xn))
            izn = argmin(abs.(z .- zn))
            @views F .= F ./ maximum(F[ixn,izn,:])
        end
    end

    isnothing(xlims[1]) ? xmin=x[1] : xmin=xlims[1]
    isnothing(xlims[2]) ? xmax=x[end] : xmax=xlims[2]
    isnothing(zlims[1]) ? zmin=z[1] : zmin=zlims[1]
    isnothing(zlims[2]) ? zmax=z[end] : zmax=zlims[2]

    fig = mak.Figure(resolution=(950,992), fontsize=14)
    if new_window
        mak.display(mak.Screen(), fig)
    end
    ax = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="z ($szu)", aspect)
    mak.xlims!(ax, (xmin,xmax))
    mak.ylims!(ax, (zmin,zmax))
    mak.display(fig)

    hm = mak.heatmap!(ax, x, z, F; colormap=cmap, colorrange=(vmin,vmax))
    mak.Colorbar(fig[2,1], hm; vertical=false, label="Time averaged |S|")

    if save
        ext = splitext(fname)[end]
        fname_fig = replace(fname, ext => ".png")
        mak.save(fname_fig, fig)
    end
    return nothing
end


# ******************************************************************************************
function poynting(Hy, Ex, Ez)
    return @. sqrt((-Ez * Hy)^2 + (Ex * Hy)^2)
end
