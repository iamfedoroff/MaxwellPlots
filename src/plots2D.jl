function plot2D(
    fname, var, t0;
    xu=nothing, zu=nothing, tu=nothing, norm=true, colormap=nothing, colorrange=nothing,
    aspect=1, new_window=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "fields/t")
    F = HDF5.read(fp, "fields/" * string(var))
    HDF5.close(fp)

    it = argmin(abs.(t .- t0))

    isnothing(xu) ? xu = space_units(x) : nothing
    isnothing(zu) ? zu = space_units(z) : nothing
    isnothing(tu) ? tu = time_units(t) : nothing
    xu = zu = max(xu, zu)
    sxu = space_units_name(xu)
    szu = space_units_name(zu)
    stu = time_units_name(tu)

    @. x = x / xu
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
    xu=nothing, zu=nothing, tu=nothing, xlims=(nothing,nothing), zlims=(nothing,nothing),
    norm=true, colormap=nothing, colorrange=nothing, aspect=1, movie=false,
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

    isnothing(xu) ? xu = space_units(x) : nothing
    isnothing(zu) ? zu = space_units(z) : nothing
    isnothing(tu) ? tu = time_units(t) : nothing
    xu = zu = max(xu, zu)
    sxu = space_units_name(xu)
    szu = space_units_name(zu)
    stu = time_units_name(tu)

    @. x = x / xu
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

    isnothing(xlims[1]) ? xmin=x[1] : xmin=xlims[1]
    isnothing(xlims[2]) ? xmax=x[end] : xmax=xlims[2]
    isnothing(zlims[1]) ? zmin=z[1] : zmin=zlims[1]
    isnothing(zlims[2]) ? zmax=z[end] : zmax=zlims[2]

    fig = mak.Figure(size=(950,992))
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
    fname, var, x0, z0; xu=nothing, zu=nothing, tu=nothing, vmin=-1, vmax=1, norm=true,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "fields/t")
    F = HDF5.read(fp, "fields/" * string(var))
    HDF5.close(fp)

    ix0 = argmin(abs.(x .- x0))
    iz0 = argmin(abs.(z .- z0))

    isnothing(xu) ? xu = space_units(x) : nothing
    isnothing(zu) ? zu = space_units(z) : nothing
    isnothing(tu) ? tu = time_units(t) : nothing
    xu = zu = max(xu, zu)
    sxu = space_units_name(xu)
    szu = space_units_name(zu)
    stu = time_units_name(tu)

    @. x = x / xu
    @. z = z / zu
    @. t = t / tu

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    fig = mak.Figure(size=(950,992))
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


function plot2D(
    fname, var;
    xu=nothing, zu=nothing, xlims=(nothing,nothing), zlims=(nothing,nothing), norm=true,
    colormap=nothing, colorrange=nothing, aspect=1, save=false, save_fname="image.png",
    new_window=false, colorbar=true,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    z = HDF5.read(fp, "z")
    if string(var) in ("Sa", "rho_end")
        F = HDF5.read(fp, string(var))
    else
        error("Unknown variable $(string(var))")
    end
    HDF5.close(fp)

    isnothing(xu) ? xu = space_units(x) : nothing
    isnothing(zu) ? zu = space_units(z) : nothing
    xu = zu = max(xu, zu)
    sxu = space_units_name(xu)
    szu = space_units_name(zu)

    @. x = x / xu
    @. z = z / zu

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    isnothing(colorrange) ? colorrange = extrema(F) : nothing
    vmin, vmax = colorrange
    isnothing(vmin) ? vmin = minimum(F) : nothing
    isnothing(vmax) ? vmax = maximum(F) : nothing
    if vmin * vmax < 0   # diverging colormap
        isnothing(colormap) ? colormap = CMAPDIV : nothing
    else
        isnothing(colormap) ? colormap = CMAP : nothing
    end

    isnothing(xlims[1]) ? xmin=x[1] : xmin=xlims[1]
    isnothing(xlims[2]) ? xmax=x[end] : xmax=xlims[2]
    isnothing(zlims[1]) ? zmin=z[1] : zmin=zlims[1]
    isnothing(zlims[2]) ? zmax=z[end] : zmax=zlims[2]

    fig = mak.Figure(size=(950,992))
    if new_window
        mak.display(mak.Screen(), fig)
    end
    ax = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="z ($szu)", aspect)
    mak.xlims!(ax, (xmin,xmax))
    mak.ylims!(ax, (zmin,zmax))
    mak.display(fig)

    hm = mak.heatmap!(ax, x, z, F; colormap, colorrange)

    if colorbar
        mak.Colorbar(fig[2,1], hm; vertical=false, label=string(var), flipaxis=false)
    end

    if save
        mak.save(save_fname, fig)
    end
    return nothing
end


# ******************************************************************************************
function poynting(Hy, Ex, Ez)
    return @. sqrt((-Ez * Hy)^2 + (Ex * Hy)^2)
end
