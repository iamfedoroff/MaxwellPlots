function inspect3D(
    fname, svar;
    xu=1, yu=1, zu=1, tu=1, norm=false, vmin=-1, vmax=1, aspect=:data,
    xlims=(nothing,nothing), ylims=(nothing,nothing), zlims=(nothing,nothing),
    cmap=:seismic, new_window=false, movie=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "t")
    F = HDF5.read(fp, svar)
    HDF5.close(fp)

    @show extrema(F)

    inspect(
        x, y, z, t, F;
        xu, yu, zu, tu, norm, vmin, vmax, aspect, xlims, ylims, zlims, cmap, new_window,
        movie,
    )
    return nothing
end


function inspect3D_xsec(
    fname, svar, x0, y0, z0;
    xu=1, yu=1, zu=1, tu=1, norm=false, norm_point=nothing, vmin=-1, vmax=1, aspect=(1,1,1),
    xlims=(nothing,nothing), ylims=(nothing,nothing), zlims=(nothing,nothing),
    cmap=:seismic, new_window=false, movie=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "t")
    F = HDF5.read(fp, svar)
    HDF5.close(fp)

    @show extrema(F)

    inspect_xsec(
        x, y, z, t, F, x0, y0, z0;
        xu, yu, zu, tu, norm, norm_point, vmin, vmax, aspect, xlims, ylims, zlims, cmap,
        new_window, movie,
    )
    return nothing
end


function inspect3D_poynting(
    fname;
    xu=1, yu=1, zu=1, tu=1, norm=false, vmin=0, vmax=1, aspect=:data,
    xlims=(nothing,nothing), ylims=(nothing,nothing), zlims=(nothing,nothing),
    cmap=mak.Reverse(:Hiroshige), new_window=false, movie=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "t")
    Hx = HDF5.read(fp, "Hx")
    Hy = HDF5.read(fp, "Hy")
    Hz = HDF5.read(fp, "Hz")
    Ex = HDF5.read(fp, "Ex")
    Ey = HDF5.read(fp, "Ey")
    Ez = HDF5.read(fp, "Ez")
    HDF5.close(fp)

    F = @. sqrt((Ey*Hz - Ez*Hy)^2 + (Ez*Hx - Ex*Hz)^2 + (Ex*Hy - Ey*Hx)^2)

    @show extrema(F)

    inspect(
        x, y, z, t, F;
        xu, yu, zu, tu, norm, vmin, vmax, aspect, xlims, ylims, zlims, cmap, new_window,
        movie,
    )
    return nothing
end


function inspect3D_poynting_xsec(
    fname, x0, y0, z0;
    xu=1, yu=1, zu=1, tu=1, norm=false, norm_point=nothing, vmin=0, vmax=1, aspect=(1,1,1),
    xlims=(nothing,nothing), ylims=(nothing,nothing), zlims=(nothing,nothing),
    cmap=mak.Reverse(:Hiroshige), new_window=false, movie=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "t")
    Hx = HDF5.read(fp, "Hx")
    Hy = HDF5.read(fp, "Hy")
    Hz = HDF5.read(fp, "Hz")
    Ex = HDF5.read(fp, "Ex")
    Ey = HDF5.read(fp, "Ey")
    Ez = HDF5.read(fp, "Ez")
    HDF5.close(fp)

    F = @. sqrt((Ey*Hz - Ez*Hy)^2 + (Ez*Hx - Ex*Hz)^2 + (Ex*Hy - Ey*Hx)^2)

    @show extrema(F)

    inspect_xsec(
        x, y, z, t, F, x0, y0, z0;
        xu, yu, zu, tu, norm, norm_point, vmin, vmax, aspect, xlims, ylims, zlims, cmap,
        new_window, movie,
    )
    return nothing
end


function plot3D_poynting_averaged(
    fname;
    xu=1, yu=1, zu=1, norm=false, vmin=0, vmax=1, aspect=:data,
    xlims=(nothing,nothing), ylims=(nothing,nothing), zlims=(nothing,nothing),
    cmap=mak.Reverse(:Hiroshige), new_window=false, save=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    F = HDF5.read(fp, "Sa")
    HDF5.close(fp)

    @show extrema(F)

    mplot(
        x, y, z, F;
        xu, yu, zu, norm, vmin, vmax, aspect, xlims, ylims, zlims, cmap, new_window, save,
    )
    return nothing
end


function plot3D_poynting_averaged_diff(
    fname1, fname2;
    xu=1, yu=1, zu=1, norm=false, vmin=-1, vmax=1, aspect=:data,
    xlims=(nothing,nothing), ylims=(nothing,nothing), zlims=(nothing,nothing),
    cmap=:seismic, new_window=false, save=false,
)
    fp = HDF5.h5open(fname1, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    F1 = HDF5.read(fp, "Sa")
    HDF5.close(fp)

    fp = HDF5.h5open(fname2, "r")
    F2 = HDF5.read(fp, "Sa")
    HDF5.close(fp)

    F = @. F2 - F1

    @show extrema(F)

    mplot(
        x, y, z, F;
        xu, yu, zu, norm, vmin, vmax, aspect, xlims, ylims, zlims, cmap, new_window, save,
    )
    return nothing
end


function plot3D_poynting_averaged_xsec(
    fname, x0, y0, z0;
    xu=1, yu=1, zu=1, norm=false, norm_point=nothing, vmin=0, vmax=1, aspect=(1,1,1),
    xlims=(nothing,nothing), ylims=(nothing,nothing), zlims=(nothing,nothing),
    cmap=mak.Reverse(:Hiroshige), new_window=false, save=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    F = HDF5.read(fp, "Sa")
    HDF5.close(fp)

    @show extrema(F)

    mplot_xsec(
       x, y, z, F, x0, y0, z0;
        xu, yu, zu, norm, norm_point, vmin, vmax, aspect, xlims, ylims, zlims, cmap,
        new_window, save,
    )
    return nothing
end


function plot3D_poynting_averaged_xsec_diff(
    fname1, fname2, x0, y0, z0;
    xu=1, yu=1, zu=1, norm=false, norm_point=nothing, vmin=0, vmax=1,
    xlims=(nothing,nothing), ylims=(nothing,nothing), zlims=(nothing,nothing),
    aspect=(1,1,1), cmap=mak.Reverse(:Hiroshige), new_window=false, save=false,
)
    fp = HDF5.h5open(fname1, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    F1 = HDF5.read(fp, "Sa")
    HDF5.close(fp)

    fp = HDF5.h5open(fname2, "r")
    F2 = HDF5.read(fp, "Sa")
    HDF5.close(fp)

    F = @. F2 - F1

    @show extrema(F)

    mplot_xsec(
       x, y, z, F, x0, y0, z0;
        xu, yu, zu, norm, norm_point, vmin, vmax, aspect, xlims, ylims, zlims, cmap,
        new_window, save,
    )
    return nothing
end
