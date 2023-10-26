function inspect3D(
    fname, var;
    xu=1, yu=1, zu=1, tu=1,
    xlims=nothing, ylims=nothing, zlims=nothing, tlims=nothing,
    norm=true, colormap=nothing, colorrange=nothing, colorbar=true, aspect=:data,
    movie=false, movie_fname=nothing, new_window=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "fields/t")
    if var in (:Ex, :Ey, :Ez, :Hx, :Hy, :Hz)
        F = HDF5.read(fp, "fields/" * string(var))
        isnothing(colormap) ? colormap = CMAPDIV : nothing
    elseif var == :poynting
        Hx = HDF5.read(fp, "fields/Hx")
        Hy = HDF5.read(fp, "fields/Hy")
        Hz = HDF5.read(fp, "fields/Hz")
        Ex = HDF5.read(fp, "fields/Ex")
        Ey = HDF5.read(fp, "fields/Ey")
        Ez = HDF5.read(fp, "fields/Ez")
        F = poynting(Hx, Hy, Hz, Ex, Ey, Ez)
        isnothing(colormap) ? colormap = CMAP : nothing
    else
        error("Wrong input varible " * string(var))
    end
    HDF5.close(fp)

    if movie && isnothing(movie_fname)
        ext = splitext(fname)[end]
        movie_fname = replace(fname, ext => ".mp4")
    end

    inspect(
        x, y, z, t, F;
        xu, yu, zu, tu, xlims, ylims, zlims, tlims, norm, colormap, colorrange, colorbar,
        aspect, movie, movie_fname, new_window,
    )
    return nothing
end


function inspect3D_xsec(
    fname, var;
    xu=1, yu=1, zu=1, tu=1,
    xcut=nothing, ycut=nothing, zcut=nothing, tlims=nothing,
    xlims=nothing, ylims=nothing, zlims=nothing,
    norm=true, colormap=nothing, colorrange=nothing, colorbar=true, aspect=(1,1,1),
    movie=false, movie_fname=nothing, new_window=false, guidelines=true,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "fields/t")
    if var in (:Ex, :Ey, :Ez, :Hx, :Hy, :Hz)
        F = HDF5.read(fp, "fields/" * string(var))
        isnothing(colormap) ? colormap = CMAPDIV : nothing
    elseif var == :poynting
        Hx = HDF5.read(fp, "fields/Hx")
        Hy = HDF5.read(fp, "fields/Hy")
        Hz = HDF5.read(fp, "fields/Hz")
        Ex = HDF5.read(fp, "fields/Ex")
        Ey = HDF5.read(fp, "fields/Ey")
        Ez = HDF5.read(fp, "fields/Ez")
        F = poynting(Hx, Hy, Hz, Ex, Ey, Ez)
        isnothing(colormap) ? colormap = CMAP : nothing
    else
        error("Wrong input varible " * string(var))
    end
    HDF5.close(fp)

    if movie && isnothing(movie_fname)
        ext = splitext(fname)[end]
        movie_fname = replace(fname, ext => ".mp4")
    end

    inspect_xsec(
        x, y, z, t, F;
        xu, yu, zu, tu, xcut, ycut, zcut, xlims, ylims, zlims, tlims, norm, colormap,
        colorrange, colorbar, aspect, movie, movie_fname, new_window, guidelines,
    )
    return nothing
end


function inspect3D_volume(
    fname, var;
    xu=1, yu=1, zu=1,
    xcut=nothing, ycut=nothing, zcut=nothing,
    xlims=nothing, ylims=nothing, zlims=nothing,
    norm=true, colormap=nothing, colorrange=nothing, colorbar=true, aspect=:data,
    save=false, save_fname=nothing, new_window=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    if var == :Sa
        F = HDF5.read(fp, "Sa")
        isnothing(colormap) ? colormap = CMAP : nothing
    else
        error("Wrong input varible " * string(var))
    end
    HDF5.close(fp)

    if save && isnothing(save_fname)
        ext = splitext(fname)[end]
        save_fname = replace(fname, ext => ".png")
    end

    inspect_volume(
        x, y, z, F;
        xu, yu, zu, xcut, ycut, zcut, xlims, ylims, zlims, norm, colormap, colorrange,
        colorbar, aspect, save, save_fname, new_window,
    )
    return nothing
end


function plot3D(
    fname, var;
    xu=1, yu=1, zu=1,
    xlims=nothing, ylims=nothing, zlims=nothing,
    norm=true, colormap=nothing, colorrange=nothing, colorbar=true, aspect=:data,
    save=false, save_fname=nothing, new_window=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    if var == :Sa
        F = HDF5.read(fp, "Sa")
        isnothing(colormap) ? colormap = CMAP : nothing
    end
    HDF5.close(fp)

    if save && isnothing(save_fname)
        ext = splitext(fname)[end]
        save_fname = replace(fname, ext => ".png")
    end

    plot_volume(
        x, y, z, F;
        xu, yu, zu, xlims, ylims, zlims, norm, colormap, colorrange, colorbar, aspect, save,
        save_fname, new_window,
    )
    return nothing
end


function plot3D_diff(
    fname1, fname2, var;
    xu=1, yu=1, zu=1,
    xlims=nothing, ylims=nothing, zlims=nothing,
    norm=true, colormap=CMAPDIV, colorrange=nothing, colorbar=true, aspect=:data,
    save=false, save_fname=nothing, new_window=false,
)
    fp = HDF5.h5open(fname1, "r")
    x1 = HDF5.read(fp, "x")
    y1 = HDF5.read(fp, "y")
    z1 = HDF5.read(fp, "z")
    if var == :Sa
        F1 = HDF5.read(fp, "Sa")
    end
    HDF5.close(fp)

    fp = HDF5.h5open(fname2, "r")
    x2 = HDF5.read(fp, "x")
    y2 = HDF5.read(fp, "y")
    z2 = HDF5.read(fp, "z")
    if var == :Sa
        F2 = HDF5.read(fp, "Sa")
    end
    HDF5.close(fp)

    @assert x1==x2 && y1==y2 && z1==z2

    F = @. F2 - F1

    if save && isnothing(save_fname)
        ext = splitext(fname2)[end]
        save_fname = replace(fname2, ext => "_diff.png")
    end

    plot_volume(
        x1, y1, z1, F;
        xu, yu, zu, xlims, ylims, zlims, norm, colormap, colorrange, colorbar, aspect, save,
        save_fname, new_window,
    )
    return nothing
end


function plot3D_xsec(
    fname, var;
    xu=1, yu=1, zu=1,
    xcut=nothing, ycut=nothing, zcut=nothing,
    xlims=nothing, ylims=nothing, zlims=nothing,
    norm=true, colormap=nothing, colorrange=nothing, colorbar=true, aspect=(1,1,1),
    save=false, save_fname=nothing, new_window=false, guidelines=true,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    if var == :Sa
        F = HDF5.read(fp, "Sa")
        isnothing(colormap) ? colormap = CMAP : nothing
    end
    HDF5.close(fp)

    if save && isnothing(save_fname)
        ext = splitext(fname)[end]
        save_fname = replace(fname, ext => ".png")
    end

    plot_volume_xsec(
        x, y, z, F;
        xu, yu, zu, xcut, ycut, zcut, xlims, ylims, zlims, norm, colormap, colorrange,
        colorbar, aspect, save, save_fname, new_window, guidelines,
    )
    return nothing
end


function plot3D_xsec_diff(
    fname1, fname2, var;
    xu=1, yu=1, zu=1,
    xcut=nothing, ycut=nothing, zcut=nothing,
    xlims=nothing, ylims=nothing, zlims=nothing,
    norm=true, colormap=CMAPDIV, colorrange=nothing, colorbar=true, aspect=(1,1,1),
    save=false, save_fname=nothing, new_window=false, guidelines=true,
)
    fp = HDF5.h5open(fname1, "r")
    x1 = HDF5.read(fp, "x")
    y1 = HDF5.read(fp, "y")
    z1 = HDF5.read(fp, "z")
    if var == :Sa
        F1 = HDF5.read(fp, "Sa")
    end
    HDF5.close(fp)

    fp = HDF5.h5open(fname2, "r")
    x2 = HDF5.read(fp, "x")
    y2 = HDF5.read(fp, "y")
    z2 = HDF5.read(fp, "z")
    if var == :Sa
        F2 = HDF5.read(fp, "Sa")
    end
    HDF5.close(fp)

    @assert x1==x2 && y1==y2 && z1==z2

    F = @. F2 - F1

    if save && isnothing(save_fname)
        ext = splitext(fname2)[end]
        save_fname = replace(fname2, ext => "_diff.png")
    end

    plot_volume_xsec(
        x1, y1, z1, F;
        xu, yu, zu, xcut, ycut, zcut, xlims, ylims, zlims, norm, colormap, colorrange,
        colorbar, aspect, save, save_fname, new_window, guidelines,
    )
    return nothing
end


# ******************************************************************************************
@kernel function poynting_kernel!(S, Hx, Hy, Hz, Ex, Ey, Ez)
    I = @index(Global)
    @inbounds begin
        Sx = Ey[I] * Hz[I] - Ez[I] * Hy[I]
        Sy = Ez[I] * Hx[I] - Ex[I] * Hz[I]
        Sz = Ex[I] * Hy[I] - Ey[I] * Hx[I]
        S[I] = sqrt(Sx^2 + Sy^2 + Sz^2)
    end
end
function poynting(Hx, Hy, Hz, Ex, Ey, Ez)
    S = similar(Hx)
    backend = get_backend(S)
    ndrange = size(S)
    poynting_kernel!(backend)(S, Hx, Hy, Hz, Ex, Ey, Ez; ndrange)
    return S
end
