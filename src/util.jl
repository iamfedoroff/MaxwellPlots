const CMAP = mak.Reverse(:Hiroshige)
const CMAPDIV = :bwr


halfint(N) = iseven(N) ? div(N,2) : div(N,2)+1


function space_units(x::AbstractArray)
    return space_units(max(extrema(x)...))
end


function space_units(x)
    if abs(x) < 1e-6
        xu = 1e-9
    elseif abs(x) < 1e-3
        xu = 1e-6
    elseif abs(x) < 1
        xu = 1e-3
    else
        xu = 1.0
    end
    return xu
end


function space_units_name(xu)
    sxu = "arb. u."
    if xu == 1
        sxu = "m"
    elseif xu == 1e-2
        sxu = "cm"
    elseif  xu == 1e-3
        sxu = "mm"
    elseif xu == 1e-6
        sxu = "μm"
    elseif xu == 1e-9
        sxu = "nm"
    end
    return sxu
end


function time_units(t::AbstractArray)
    return time_units(max(extrema(t)...))
end


function time_units(t)
    if abs(t) < 1e-12
        tu = 1e-15
    elseif abs(t) < 1e-9
        tu = 1e-12
    elseif abs(t) < 1e-6
        tu = 1e-9
    elseif abs(t) < 1e-3
        tu = 1e-6
    elseif abs(t) < 1
        tu = 1e-3
    else
        tu = 1.0
    end
    return tu
end


function time_units_name(tu)
    stu = "arb. u."
    if tu == 1
        stu = "s"
    elseif  tu == 1e-3
        stu = "ms"
    elseif tu == 1e-6
        stu = "μs"
    elseif tu == 1e-9
        stu = "ns"
    elseif tu == 1e-12
        stu = "ps"
    elseif tu == 1e-15
        stu = "fs"
    elseif tu == 1e-18
        stu = "as"
    end
    return stu
end


function apply_limits(x, y, z, F; xlims=nothing, ylims=nothing, zlims=nothing)
    isnothing(xlims) ? xlims=(nothing,nothing) : nothing
    isnothing(ylims) ? ylims=(nothing,nothing) : nothing
    isnothing(zlims) ? zlims=(nothing,nothing) : nothing
    xmin, xmax = xlims
    ymin, ymax = ylims
    zmin, zmax = zlims
    isnothing(xmin) ? xmin=x[1] : nothing
    isnothing(xmax) ? xmax=x[end] : nothing
    isnothing(ymin) ? ymin=y[1] : nothing
    isnothing(ymax) ? ymax=y[end] : nothing
    isnothing(zmin) ? zmin=z[1] : nothing
    isnothing(zmax) ? zmax=z[end] : nothing
    ix1, ix2 = argmin(abs.(x.-xmin)), argmin(abs.(x.-xmax))
    iy1, iy2 = argmin(abs.(y.-ymin)), argmin(abs.(y.-ymax))
    iz1, iz2 = argmin(abs.(z.-zmin)), argmin(abs.(z.-zmax))
    return x[ix1:ix2], y[iy1:iy2], z[iz1:iz2], F[ix1:ix2,iy1:iy2,iz1:iz2]
end


function apply_limits(
    x, y, z, t, F; xlims=nothing, ylims=nothing, zlims=nothing, tlims=nothing,
)
    isnothing(xlims) ? xlims=(nothing,nothing) : nothing
    isnothing(ylims) ? ylims=(nothing,nothing) : nothing
    isnothing(zlims) ? zlims=(nothing,nothing) : nothing
    isnothing(tlims) ? tlims=(nothing,nothing) : nothing
    xmin, xmax = xlims
    ymin, ymax = ylims
    zmin, zmax = zlims
    tmin, tmax = tlims
    isnothing(xmin) ? xmin=x[1] : nothing
    isnothing(xmax) ? xmax=x[end] : nothing
    isnothing(ymin) ? ymin=y[1] : nothing
    isnothing(ymax) ? ymax=y[end] : nothing
    isnothing(zmin) ? zmin=z[1] : nothing
    isnothing(zmax) ? zmax=z[end] : nothing
    isnothing(tmin) ? tmin=t[1] : nothing
    isnothing(tmax) ? tmax=t[end] : nothing
    ix1, ix2 = argmin(abs.(x.-xmin)), argmin(abs.(x.-xmax))
    iy1, iy2 = argmin(abs.(y.-ymin)), argmin(abs.(y.-ymax))
    iz1, iz2 = argmin(abs.(z.-zmin)), argmin(abs.(z.-zmax))
    it1, it2 = argmin(abs.(t.-tmin)), argmin(abs.(t.-tmax))
    return x[ix1:ix2], y[iy1:iy2], z[iz1:iz2], t[it1:it2],
           F[ix1:ix2,iy1:iy2,iz1:iz2,it1:it2]
end


function plot_geometry(
    x, z, geometry;
    xu=nothing, zu=nothing, cmap=CMAP, aspect=1,
    smooth_interfaces=false, new_window=false,
)
    if typeof(geometry) <: Function
        Nx, Nz = length(x), length(z)
        F = [geometry(x[ix],z[iz]) ? 1 : 0 for ix=1:Nx, iz=1:Nz]
    else
        F = geometry
    end

    if smooth_interfaces
        F = moving_average(F, 2)
    end

    isnothing(xu) ? xu = space_units(x) : nothing
    isnothing(zu) ? zu = space_units(z) : nothing
    xu = zu = max(xu, zu)
    sxu = space_units_name(xu)
    szu = space_units_name(zu)

    xx = x / xu
    zz = z / zu

    fig = mak.Figure(size=(950,992))
    if new_window
        mak.display(mak.Screen(), fig)
    end
    ax = mak.Axis(
        fig[1,1];
        xlabel="x ($sxu)", ylabel="z ($szu)", title="Material geometry", aspect,
    )
    mak.display(fig)
    hm = mak.heatmap!(
        ax, xx, zz, F; colormap=cmap, colorrange=(0,1),
    )
    # mak.Colorbar(fig[2,1], hm; vertical=false, label="geometry")
    return nothing
end


function plot_geometry(
    x, y, z, geometry;
    xu=nothing, yu=nothing, zu=nothing, xlims=nothing, ylims=nothing, zlims=nothing,
    colormap=CMAP, aspect=:data, smooth_interfaces=false, new_window=false, algorithm=:iso,
    isovalue=1, absorption=1, save=false, save_fname="out.png",
)
    if typeof(geometry) <: Function
        Nx, Ny, Nz = length(x), length(y), length(z)
        F = [geometry(x[ix],y[iy],z[iz]) ? 1 : 0 for ix=1:Nx, iy=1:Ny, iz=1:Nz]
    else
        F = geometry
    end

    x, y, z, F = apply_limits(x, y, z, F; xlims, ylims, zlims)

    if smooth_interfaces
        F = moving_average(F, 2)
    end

    isnothing(xu) ? xu = space_units(x) : nothing
    isnothing(yu) ? yu = space_units(y) : nothing
    isnothing(zu) ? zu = space_units(z) : nothing
    xu = yu = zu = max(xu, yu, zu)
    sxu = space_units_name(xu)
    syu = space_units_name(yu)
    szu = space_units_name(zu)

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


function plot_waveform(model; tu=1)
    (; field, sources, t) = model
    (; waveform, p, icomp) = sources[1]

    isnothing(tu) ? tu = time_units(t) : nothing
    stu = time_units_name(tu)

    tt = t / tu

    component = fieldnames(typeof(field))[icomp]

    # F = @. waveform(0, 0, 0, t, (p,))
    F = @. waveform(0, 0, 0, t)

    fig = mak.Figure(size=(950,992))
    mak.display(mak.Screen(), fig)
    ax = mak.Axis(fig[1,1]; xlabel="t ($stu)", ylabel=string(component))
    mak.lines!(ax, tt, F)
    mak.display(fig)
    return nothing
end


# https://discourse.julialang.org/t/makie-figure-resolution-makie-primary-resolution-deprecated/93854/4
function primary_resolution()
    monitor = mak.GLFW.GetPrimaryMonitor()
    videomode = mak.MonitorProperties(monitor).videomode
    return (videomode.width, videomode.height)
end
