const CMAP = mak.Reverse(:Hiroshige)
const CMAPDIV = :bwr


halfint(N) = iseven(N) ? div(N,2) : div(N,2)+1


function default_colormap(F)
    return prod(extrema(F)) < 0 ? CMAPDIV : CMAP
end


function default_colorrange(F)
    ex = extrema(F)
    if prod(ex) < 0
        exmax = max(ex...)
        colorrange = (-exmax, exmax)
    else
        colorrange = ex
    end
    return colorrange
end


function units(x::AbstractArray)
    return units(max(maximum(abs,x)))
end


function units(x)
    if abs(x) < 1e-12
        xu = 1e-15
    elseif abs(x) < 1e-9
        xu = 1e-12
    elseif abs(x) < 1e-6
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


function units_name_space(xu)
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


function units_name_time(tu)
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


function indices_of_limits(x, xlims)
    xlims = isnothing(xlims) ? (nothing,nothing) : xlims
    xmin, xmax = xlims
    xmin = isnothing(xmin) ? x[1] : xmin
    xmax = isnothing(xmax) ? x[end] : xmax
    ix1, ix2 = argmin(abs.(x.-xmin)), argmin(abs.(x.-xmax))
    return ix1, ix2
end


function apply_limits(x, y, F; xlims=nothing, ylims=nothing)
    ix1, ix2 = indices_of_limits(x, xlims)
    iy1, iy2 = indices_of_limits(y, ylims)
    return x[ix1:ix2], y[iy1:iy2], F[ix1:ix2,iy1:iy2]
end


function apply_limits(x, y, z, F; xlims=nothing, ylims=nothing, zlims=nothing)
    ix1, ix2 = indices_of_limits(x, xlims)
    iy1, iy2 = indices_of_limits(y, ylims)
    iz1, iz2 = indices_of_limits(z, zlims)
    return x[ix1:ix2], y[iy1:iy2], z[iz1:iz2], F[ix1:ix2,iy1:iy2,iz1:iz2]
end


function apply_limits(
    x, y, z, t, F; xlims=nothing, ylims=nothing, zlims=nothing, tlims=nothing,
)
    ix1, ix2 = indices_of_limits(x, xlims)
    iy1, iy2 = indices_of_limits(y, ylims)
    iz1, iz2 = indices_of_limits(z, zlims)
    it1, it2 = indices_of_limits(t, tlims)
    return x[ix1:ix2], y[iy1:iy2], z[iz1:iz2], t[it1:it2],
           F[ix1:ix2,iy1:iy2,iz1:iz2,it1:it2]
end


function grid_dims(fname)
    fp = HDF5.h5open(fname, "r")
    fpkeys = keys(fp)
    if "y" in fpkeys
        n = 3
    elseif "x" in fpkeys
        n = 2
    else
        n = 1
    end
    HDF5.close(fp)
    return n
end


function plot_waveform(model; tu=1)
    (; field, sources, t) = model
    (; waveform, p, icomp) = sources[1]

    isnothing(tu) ? tu = units(t) : nothing
    stu = units_name_time(tu)

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
