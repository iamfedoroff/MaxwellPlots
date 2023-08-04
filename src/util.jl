function space_units_string(xu)
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


function time_units_string(tu)
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


function plot_geometry(
    x, z, geometry;
    xu=1, zu=1, cmap=mak.Reverse(:Hiroshige), aspect=1,
    smooth_interfaces=false, new_window=false,
)
    Nx, Nz = length(x), length(z)
    F = [geometry(x[ix],z[iz]) ? 1 : 0 for ix=1:Nx, iz=1:Nz]

    if smooth_interfaces
        F = moving_average(F, 2)
    end

    xx = x / xu
    zz = z / zu
    sxu = space_units_string(xu)
    szu = space_units_string(zu)

    fig = mak.Figure(resolution=(950,992), fontsize=14)
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
    xu=1, yu=1, zu=1, cmap=mak.Reverse(:Hiroshige), aspect=:data,
    smooth_interfaces=false, new_window=false,
)
    Nx, Ny, Nz = length(x), length(y), length(z)
    F = [geometry(x[ix],y[iy],z[iz]) ? 1 : 0 for ix=1:Nx, iy=1:Ny, iz=1:Nz]

    if smooth_interfaces
        F = moving_average(F, 2)
    end

    xx = x / xu
    yy = y / yu
    zz = z / zu
    sxu = space_units_string(xu)
    syu = space_units_string(yu)
    szu = space_units_string(zu)

    fig = mak.Figure(resolution=(950,992), fontsize=14)
    if new_window
        mak.display(mak.Screen(), fig)
    end
    ax = mak.Axis3(
        fig[1,1];
        xlabel="x ($sxu)", ylabel="y ($syu)", zlabel="z ($szu)",
        title="Material geometry", aspect, perspectiveness=0,
    )
    img = mak.volume!(
        ax, xx, yy, zz, F; colormap=cmap, algorithm=:iso, isovalue=1,
    )
    # mak.Colorbar(fig[1,2], img; label="geometry")
    mak.display(fig)
    return nothing
end


function plot_waveform(model; tu=1)
    (; field, source, t) = model
    (; waveform, p, icomp) = source

    tt = t / tu
    stu = time_units_string(tu)

    component = fieldnames(typeof(field))[icomp]

    F = @. waveform(t, (p,))

    fig = mak.Figure(resolution=(950,992), fontsize=14)
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
