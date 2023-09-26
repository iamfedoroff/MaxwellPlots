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
    t = HDF5.read(fp, "fields/t")
    F = HDF5.read(fp, "fields/"*svar)
    HDF5.close(fp)

    @show extrema(F)

    ext = splitext(fname)[end]
    movie_fname = replace(fname, ext => ".mp4")

    inspect(
        x, y, z, t, F;
        xu, yu, zu, tu, norm, vmin, vmax, aspect, xlims, ylims, zlims, cmap, new_window,
        movie, movie_fname,
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
    t = HDF5.read(fp, "fields/t")
    F = HDF5.read(fp, "fields/"*svar)
    HDF5.close(fp)

    @show extrema(F)

    ext = splitext(fname)[end]
    movie_fname = replace(fname, ext => ".mp4")

    inspect_xsec(
        x, y, z, t, F, x0, y0, z0;
        xu, yu, zu, tu, norm, norm_point, vmin, vmax, aspect, xlims, ylims, zlims, cmap,
        new_window, movie, movie_fname,
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
    t = HDF5.read(fp, "fields/t")
    Hx = HDF5.read(fp, "fields/Hx")
    Hy = HDF5.read(fp, "fields/Hy")
    Hz = HDF5.read(fp, "fields/Hz")
    Ex = HDF5.read(fp, "fields/Ex")
    Ey = HDF5.read(fp, "fields/Ey")
    Ez = HDF5.read(fp, "fields/Ez")
    HDF5.close(fp)

    F = poynting(Hx, Hy, Hz, Ex, Ey, Ez)
    @show extrema(F)

    ext = splitext(fname)[end]
    movie_fname = replace(fname, ext => ".mp4")

    inspect(
        x, y, z, t, F;
        xu, yu, zu, tu, norm, vmin, vmax, aspect, xlims, ylims, zlims, cmap, new_window,
        movie, movie_fname
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
    t = HDF5.read(fp, "fields/t")
    Hx = HDF5.read(fp, "fields/Hx")
    Hy = HDF5.read(fp, "fields/Hy")
    Hz = HDF5.read(fp, "fields/Hz")
    Ex = HDF5.read(fp, "fields/Ex")
    Ey = HDF5.read(fp, "fields/Ey")
    Ez = HDF5.read(fp, "fields/Ez")
    HDF5.close(fp)

    F = poynting(Hx, Hy, Hz, Ex, Ey, Ez)
    @show extrema(F)

    ext = splitext(fname)[end]
    movie_fname = replace(fname, ext => ".mp4")

    inspect_xsec(
        x, y, z, t, F, x0, y0, z0;
        xu, yu, zu, tu, norm, norm_point, vmin, vmax, aspect, xlims, ylims, zlims, cmap,
        new_window, movie, movie_fname,
    )
    return nothing
end


function inspect3D_poynting_averaged(
    fname; xu=1, yu=1, zu=1, colormap=mak.Reverse(:Hiroshige), colorrange=(0,1),
    new_window=false, colorbar=true,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    F = HDF5.read(fp, "Sa")
    HDF5.close(fp)

    Nx, Ny, Nz = size(F)
    ix = iseven(Nx) ? div(Nx,2) : div(Nx,2)+1
    iy = iseven(Ny) ? div(Ny,2) : div(Ny,2)+1
    iz = iseven(Nz) ? div(Nz,2) : div(Nz,2)+1
    Lx = x[end] - x[1]
    Ly = y[end] - y[1]
    Lz = z[end] - z[1]

    @. x = x / xu
    @. y = y / yu
    @. z = z / zu

    @show extrema(F)
    F .= F ./ maximum(F)

    # --------------------------------------------------------------------------------------
    fig = mak.Figure(resolution=(950,992), fontsize=14)

    ax = mak.Axis3(
        fig[1,1]; perspectiveness=0, xlabel="x (μm)", ylabel="y (μm)", zlabel="z (μm)",
        aspect=(Lx/Lz,Ly/Lz,1),
    )
    mak.xlims!(ax, (x[1],x[end]))
    mak.ylims!(ax, (y[1],y[end]))
    mak.zlims!(ax, (z[1],z[end]))
    ax.title[] = @sprintf("%.3f μm  %.3f μm  %.3f μm", x[ix], y[iy], z[iz])

    vol = mak.volumeslices!(
        ax, x, y, z, F; bbox_visible=false, colormap, colorrange,
    )
    vol.update_yz[](ix)
    vol.update_xz[](iy)
    vol.update_xy[](iz)

    if colorbar
        mak.Colorbar(fig[1,2], vol.heatmap_xy[])
    end

    # --------------------------------------------------------------------------------------
    sg = mak.SliderGrid(
        fig[2,1],
        (label="x", range=1:length(x), startvalue=ix),
        (label="y", range=1:length(y), startvalue=iy),
        (label="z", range=1:length(z), startvalue=iz),
    )
    mak.on(sg.sliders[1].value) do i
        ix = i
        vol.update_yz[](ix)
        ax.title[] = @sprintf("%.3f μm  %.3f μm  %.3f μm", x[ix], y[iy], z[iz])
    end
    mak.on(sg.sliders[2].value) do i
        iy = i
        vol.update_xz[](iy)
        ax.title[] = @sprintf("%.3f μm  %.3f μm  %.3f μm", x[ix], y[iy], z[iz])
    end
    mak.on(sg.sliders[3].value) do i
        iz = i
        vol.update_xy[](iz)
        ax.title[] = @sprintf("%.3f μm  %.3f μm  %.3f μm", x[ix], y[iy], z[iz])
    end

    # hmaps = [vol.heatmap_yz[], vol.heatmap_xz[], vol.heatmap_xy[]]
    # toggles = [mak.Toggle(sg.layout[i,4], active=true) for i in 1:length(hmaps)]
    # map(zip(hmaps, toggles)) do (h, t)
    #     mak.connect!(h.visible, t.active)
    # end

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end

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
    cmap=mak.Reverse(:Hiroshige), new_window=false, save=false, guidelines=true,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    F = HDF5.read(fp, "Sa")
    HDF5.close(fp)

    @show extrema(F)

    ext = splitext(fname)[end]
    fname_fig = replace(fname, ext => ".png")

    mplot_xsec(
       x, y, z, F, x0, y0, z0;
       xu, yu, zu, norm, norm_point, vmin, vmax, aspect, xlims, ylims, zlims, cmap,
       new_window, save, fname_fig, guidelines,
    )
    return nothing
end


function plot3D_poynting_averaged_xsec_diff(
    fname1, fname2, x0, y0, z0;
    xu=1, yu=1, zu=1, norm=false, norm_point=nothing, vmin=0, vmax=1,
    xlims=(nothing,nothing), ylims=(nothing,nothing), zlims=(nothing,nothing),
    aspect=(1,1,1), cmap=mak.Reverse(:Hiroshige), new_window=false, save=false,
    guidelines=true,
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

    ext = splitext(fname1)[end]
    fname_fig = replace(fname1, ext => "_diff.png")

    mplot_xsec(
       x, y, z, F, x0, y0, z0;
        xu, yu, zu, norm, norm_point, vmin, vmax, aspect, xlims, ylims, zlims, cmap,
        new_window, save, fname_fig, guidelines,
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
