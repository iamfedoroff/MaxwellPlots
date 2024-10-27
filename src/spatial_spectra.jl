function spectral_grid(x, y)
    Nx = length(x)
    Ny = length(y)
    dx = x[2] - x[1]
    dy = y[2] - y[1]
    kx = 2*pi * FFTW.fftfreq(Nx, 1/dx)
    ky = 2*pi * FFTW.fftfreq(Ny, 1/dy)
    kx = FFTW.ifftshift(kx)
    ky = FFTW.ifftshift(ky)
    isodd(Nx) ? kx = kx[2:end] : nothing
    isodd(Ny) ? ky = ky[2:end] : nothing
    return kx, ky
end


function spectrum(F; func=abs2)
    S = FFTW.fft(F)
    S = @. func(S)
    S[1,1] = eps(eltype(S))   # supress DC component
    S = FFTW.ifftshift(S)
    Nx, Ny = size(F)
    isodd(Nx) ? S = S[2:end,:] : nothing
    isodd(Ny) ? S = S[:,2:end] : nothing
    return S
end


function plot_spatial_spectrum(
    fname, var, zcut;
    gmask = true,
    nopml = true,
    xu = nothing,
    yu = nothing,
    zu = nothing,
    kxu = nothing,
    kyu = nothing,
    xmin = nothing,
    xmax = nothing,
    ymin = nothing,
    ymax = nothing,
    kxmin = nothing,
    kxmax = nothing,
    kymin = nothing,
    kymax = nothing,
    colormapF = :linear_worb_100_25_c53_n256,
    colormapS = :dense,
    colorrangeF = (0,1),
    colorrangeS = (0,1),
    colorscaleS = identity,
    save = false,
    save_fname = "out.png",
    new_window = false,
)
    fp = HDF5.h5open(fname, "r")
    if nopml
        ix1, ix2, iy1, iy2, iz1, iz2 = HDF5.read(fp, "pml")
        x = fp["x"][ix1:ix2]
        y = fp["y"][iy1:iy2]
        z = fp["z"][iz1:iz2]
        if gmask
            G = fp["geometry"][ix1:ix2,iy1:iy2,iz1:iz2]
        end
        F = fp[string(var)][ix1:ix2,iy1:iy2,iz1:iz2]
    else
        x = HDF5.read(fp, "x")
        y = HDF5.read(fp, "y")
        z = HDF5.read(fp, "z")
        if gmask
            G = HDF5.read(fp, "geometry")
        end
        F = HDF5.read(fp, string(var))
    end
    HDF5.close(fp)

    xu = isnothing(xu) ? units(x) : xu
    yu = isnothing(yu) ? units(y) : yu
    zu = isnothing(zu) ? units(z) : zu
    kxu = isnothing(kxu) ? 1/xu : kxu
    kyu = isnothing(kyu) ? 1/yu : kyu
    sxu = units_name_space(xu)
    syu = units_name_space(yu)
    szu = units_name_space(zu)

    # --------------------------------------------------------------------------------------
    if gmask
        @. F = F * G
    end

    iz0 = argmin(abs.(z .- zcut))
    Fxy = F[:,:,iz0]
    Fxy .= Fxy ./ maximum(Fxy)   # normalize

    kx, ky = spectral_grid(x, y)
    Sxy = spectrum(Fxy)
    Sxy .= Sxy ./ maximum(Sxy)   # normalize

    # --------------------------------------------------------------------------------------
    fig = mak.Figure(size=(1400,800))

    ax1 = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="y ($syu)", title="z = $(z[iz0]/zu) ($szu)")
    hm = mak.heatmap!(ax1, x/xu, y/yu, Fxy; colormap=colormapF, colorrange=colorrangeF)
    mak.Colorbar(fig[2,1], hm; vertical=false, flipaxis=false, label="$(string(var))")

    ax2 = mak.Axis(fig[1,2]; xlabel="kx", ylabel="ky")
    hm = mak.heatmap!(ax2, kx/kxu, ky/kyu, Sxy; colormap=colormapS, colorrange=colorrangeS, colorscale=colorscaleS)
    mak.Colorbar(fig[2,2], hm; vertical=false, flipaxis=false, label="Power spectrum (arb.u.)")

    xmin = isnothing(xmin) ? x[1] : xmin
    xmax = isnothing(xmax) ? x[end] : xmax
    ymin = isnothing(ymin) ? y[1] : ymin
    ymax = isnothing(ymax) ? y[end] : ymax
    kxmin = isnothing(kxmin) ? kx[1] : kxmin
    kxmax = isnothing(kxmax) ? kx[end] : kxmax
    kymin = isnothing(kymin) ? ky[1] : kymin
    kymax = isnothing(kymax) ? ky[end] : kymax
    mak.xlims!(ax1, (xmin/xu, xmax/xu))
    mak.ylims!(ax1, (ymin/xu, ymax/xu))
    mak.xlims!(ax2, (kxmin/kxu, kxmax/kxu))
    mak.ylims!(ax2, (kymin/kxu, kymax/kxu))

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


function inspect_spatial_spectrum(
    fname, var;
    zcut = nothing,
    gmask = true,
    nopml = true,
    xu = nothing,
    yu = nothing,
    zu = nothing,
    kxu = nothing,
    kyu = nothing,
    xmin = nothing,
    xmax = nothing,
    ymin = nothing,
    ymax = nothing,
    kxmin = nothing,
    kxmax = nothing,
    kymin = nothing,
    kymax = nothing,
    colormapF = :linear_worb_100_25_c53_n256,
    colormapS = :dense,
    colorrangeF = (0,1),
    colorrangeS = (0,1),
    colorscaleS = identity,
    save = false,
    save_fname = "out.png",
    new_window = false,
)
    fp = HDF5.h5open(fname, "r")
    if nopml
        ix1, ix2, iy1, iy2, iz1, iz2 = HDF5.read(fp, "pml")
        x = fp["x"][ix1:ix2]
        y = fp["y"][iy1:iy2]
        z = fp["z"][iz1:iz2]
        if gmask
            G = fp["geometry"][ix1:ix2,iy1:iy2,iz1:iz2]
        end
        F = fp[string(var)][ix1:ix2,iy1:iy2,iz1:iz2]
    else
        x = HDF5.read(fp, "x")
        y = HDF5.read(fp, "y")
        z = HDF5.read(fp, "z")
        if gmask
            G = HDF5.read(fp, "geometry")
        end
        F = HDF5.read(fp, string(var))
    end
    HDF5.close(fp)

    xu = isnothing(xu) ? units(x) : xu
    yu = isnothing(yu) ? units(y) : yu
    zu = isnothing(zu) ? units(z) : zu
    kxu = isnothing(kxu) ? 1/xu : kxu
    kyu = isnothing(kyu) ? 1/yu : kyu
    sxu = units_name_space(xu)
    syu = units_name_space(yu)
    szu = units_name_space(zu)

    # --------------------------------------------------------------------------------------
    if gmask
        @. F = F * G
    end

    F .= F ./ maximum(F)

    kx, ky = spectral_grid(x, y)
    S1 = spectrum(F[:,:,1])

    S = zeros(size(S1)..., length(z))
    @. S[:,:,1] = S1
    for iz=2:length(z)
        S[:,:,iz] .= spectrum(F[:,:,iz])
    end
    S .= S ./ maximum(S)

    # --------------------------------------------------------------------------------------
    zcut = isnothing(zcut) ? z[1] + (z[end] - z[1]) / 2 : zcut
    iz = argmin(abs.(z .- zcut))

    fig = mak.Figure(size=(1400,800))

    ax1 = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="y ($syu)")
    hm1 = mak.heatmap!(ax1, x/xu, y/yu, F[:,:,iz]; colormap=colormapF, colorrange=colorrangeF)
    mak.Colorbar(fig[2,1], hm1; vertical=false, flipaxis=false, label="$(string(var))")
    ax1.title[] = @sprintf("%d:     %.3f (%s)", iz, z[iz]/zu, szu)

    ax2 = mak.Axis(fig[1,2]; xlabel="kx", ylabel="ky")
    hm2 = mak.heatmap!(ax2, kx/kxu, ky/kyu, S[:,:,iz]; colormap=colormapS, colorrange=colorrangeS, colorscale=colorscaleS)
    mak.Colorbar(fig[2,2], hm2; vertical=false, flipaxis=false, label="Power spectrum (arb.u.)")

    xmin = isnothing(xmin) ? x[1] : xmin
    xmax = isnothing(xmax) ? x[end] : xmax
    ymin = isnothing(ymin) ? y[1] : ymin
    ymax = isnothing(ymax) ? y[end] : ymax
    kxmin = isnothing(kxmin) ? kx[1] : kxmin
    kxmax = isnothing(kxmax) ? kx[end] : kxmax
    kymin = isnothing(kymin) ? ky[1] : kymin
    kymax = isnothing(kymax) ? ky[end] : kymax
    mak.xlims!(ax1, (xmin/xu, xmax/xu))
    mak.ylims!(ax1, (ymin/xu, ymax/xu))
    mak.xlims!(ax2, (kxmin/kxu, kxmax/kxu))
    mak.ylims!(ax2, (kymin/kxu, kymax/kxu))

    sg = mak.SliderGrid(fig[3,1:2], (label="z", range=1:length(z), startvalue=iz))
    mak.on(sg.sliders[1].value) do iz
        hm1[3] = F[:,:,iz]
        hm2[3] = S[:,:,iz]
        ax1.title[] = @sprintf("%d:     %.3f (%s)", iz, z[iz]/zu, szu)
    end

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


function plot_integrated_spatial_spectrum(
    fname, var;
    gmask = true,
    nopml = true,
    xu = nothing,
    yu = nothing,
    kxu = nothing,
    kyu = nothing,
    xmin = nothing,
    xmax = nothing,
    ymin = nothing,
    ymax = nothing,
    kxmin = nothing,
    kxmax = nothing,
    kymin = nothing,
    kymax = nothing,
    colormapF = :linear_worb_100_25_c53_n256,
    colormapS = :dense,
    colorrangeF = (0,1),
    colorrangeS = (0,1),
    colorscaleS = identity,
    save = false,
    save_fname = "out.png",
    new_window = false,
)
    fp = HDF5.h5open(fname, "r")
    if nopml
        ix1, ix2, iy1, iy2, iz1, iz2 = HDF5.read(fp, "pml")
        x = fp["x"][ix1:ix2]
        y = fp["y"][iy1:iy2]
        z = fp["z"][iz1:iz2]
        if gmask
            G = fp["geometry"][ix1:ix2,iy1:iy2,iz1:iz2]
        end
        F = fp[string(var)][ix1:ix2,iy1:iy2,iz1:iz2]
    else
        x = HDF5.read(fp, "x")
        y = HDF5.read(fp, "y")
        z = HDF5.read(fp, "z")
        if gmask
            G = HDF5.read(fp, "geometry")
        end
        F = HDF5.read(fp, string(var))
    end
    HDF5.close(fp)

    xu = isnothing(xu) ? units(x) : xu
    yu = isnothing(yu) ? units(y) : yu
    kxu = isnothing(kxu) ? 1/xu : kxu
    kyu = isnothing(kyu) ? 1/yu : kyu
    sxu = units_name_space(xu)
    syu = units_name_space(yu)

    # --------------------------------------------------------------------------------------
    if gmask
        @. F = F * G
    end

    F = dropdims(sum(F; dims=3); dims=3)
    F .= F ./ maximum(F)   # normalize

    kx, ky = spectral_grid(x, y)
    S = spectrum(F)
    S .= S ./ maximum(S)   # normalize

    # --------------------------------------------------------------------------------------
    fig = mak.Figure(size=(1400,800))

    ax1 = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="y ($syu)")
    hm = mak.heatmap!(ax1, x/xu, y/yu, F; colormap=colormapF, colorrange=colorrangeF)
    mak.Colorbar(fig[2,1], hm; vertical=false, flipaxis=false, label="Integrated $(string(var))")

    ax2 = mak.Axis(fig[1,2]; xlabel="kx", ylabel="ky")
    hm = mak.heatmap!(ax2, kx/kxu, ky/kyu, S; colormap=colormapS, colorrange=colorrangeS, colorscale=colorscaleS)
    mak.Colorbar(fig[2,2], hm; vertical=false, flipaxis=false, label="Integrated power spectrum (arb.u.)")

    xmin = isnothing(xmin) ? x[1] : xmin
    xmax = isnothing(xmax) ? x[end] : xmax
    ymin = isnothing(ymin) ? y[1] : ymin
    ymax = isnothing(ymax) ? y[end] : ymax
    kxmin = isnothing(kxmin) ? kx[1] : kxmin
    kxmax = isnothing(kxmax) ? kx[end] : kxmax
    kymin = isnothing(kymin) ? ky[1] : kymin
    kymax = isnothing(kymax) ? ky[end] : kymax
    mak.xlims!(ax1, (xmin/xu, xmax/xu))
    mak.ylims!(ax1, (ymin/xu, ymax/xu))
    mak.xlims!(ax2, (kxmin/kxu, kxmax/kxu))
    mak.ylims!(ax2, (kymin/kxu, kymax/kxu))

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
