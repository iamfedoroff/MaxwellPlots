function inspect3D(
    fname, svar;
    xu=1, yu=1, zu=1, tu=1, vmin=-1, vmax=1, norm=false,
    xlims=(nothing,nothing), ylims=(nothing,nothing), zlims=(nothing,nothing),
    cmap=:seismic, aspect=(1,1,1),
    movie=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "t")
    F = HDF5.read(fp, svar)
    HDF5.close(fp)

    @. x = x / xu
    @. y = y / yu
    @. z = z / zu
    @. t = t / tu
    sxu = space_units_string(xu)
    syu = space_units_string(yu)
    szu = space_units_string(zu)
    stu = time_units_string(tu)

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    # for iy=1:Ny, ix=1:Nx
    #     # if (x[ix] < 0) && (y[iy] < 0)
    #     if x[ix] < 0
    #         @. F[ix,iy,:,:] = 0
    #     end
    # end

    isnothing(xlims[1]) ? xmin=x[1] : xmin=xlims[1]
    isnothing(xlims[2]) ? xmax=x[end] : xmax=xlims[2]
    isnothing(ylims[1]) ? ymin=y[1] : ymin=ylims[1]
    isnothing(ylims[2]) ? ymax=y[end] : ymax=ylims[2]
    isnothing(zlims[1]) ? zmin=z[1] : zmin=zlims[1]
    isnothing(zlims[2]) ? zmax=z[end] : zmax=zlims[2]

    fig = mak.Figure(resolution=(950,992), fontsize=14)
    ax = mak.Axis3(
        fig[1,1];
        aspect, perspectiveness=0,
        xlabel="x ($sxu)", ylabel="y ($syu)", zlabel="z ($szu)",
    )
    mak.xlims!(ax, (xmin,xmax))
    mak.ylims!(ax, (ymin,ymax))
    mak.zlims!(ax, (zmin,zmax))
    mak.display(fig)

    cmap = mak.to_colormap(cmap)
    Nmid = div(length(cmap),2)
    # cmap[Nmid+1] = mak.RGBAf(0,0,0,0)
    @. cmap[Nmid:Nmid+2] = mak.RGBAf(0,0,0,0)

    # # alpha = [ones(Nmid-1); [0,0,0]; ones(Nmid-1)]
    # Nclip = 7
    # alpha = [ones(Nmid-Nclip); zeros(2*Nclip+1); ones(Nmid-Nclip)]
    # # alpha = [reverse(range(0.0,1.0, Nmid-1)); [0,0,0]; range(0.0,1.0, Nmid-1)]
    # cmap = [mak.RGBAf(mak.RGBf(cmap[i]), alpha[i]) for i in eachindex(cmap)]


    it = 1

    img = mak.volume!(
        ax, x, y, z, F[:,:,:,it];
        colormap=cmap, colorrange=(vmin,vmax),
        algorithm=:absorption, absorption=4f0,
        # algorithm=:iso, isovalue=0.1*vmax,
    )

    # img = mak.contour!(
    #     ax, x, y, z, F[:,:,:,it];
    #     levels=[0.1*vmin,0.1*vmax], colormap=cmap, colorrange=(vmin,vmax),
    #     alpha=1,
    # )

    mak.Colorbar(fig[1,2], img; label=svar)
    ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)

    if movie
        ext = splitext(fname)[end]
        fname_movie = replace(fname, ext => ".mp4")
        mak.record(fig, fname_movie, 1:length(t); framerate=12) do it
            img[4] = F[:,:,:,it]
            ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
            # ax.azimuth[] = 1.7pi + 0.3 * sin(2pi * it / 120)
            # ax.elevation[] = 1.7pi + 0.3 * sin(2pi * it / 120)
        end
    else
        sg = mak.SliderGrid(fig[2,1], (label="Time", range=1:length(t), startvalue=1))
        mak.on(sg.sliders[1].value) do it
            img[4] = F[:,:,:,it]
            ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end
    end
    return nothing
end


function inspect3D_xsec(
    fname, svar, x0, y0, z0;
    xu=1, yu=1, zu=1, tu=1, vmin=-1, vmax=1, norm=false,
    xlims=(nothing,nothing), ylims=(nothing,nothing), zlims=(nothing,nothing),
    cmap=:seismic,
    movie=false,
)
    fp = HDF5.h5open(fname, "r")
    x = HDF5.read(fp, "x")
    y = HDF5.read(fp, "y")
    z = HDF5.read(fp, "z")
    t = HDF5.read(fp, "t")
    F = HDF5.read(fp, svar)
    HDF5.close(fp)

    @. x = x / xu
    @. y = y / yu
    @. z = z / zu
    @. t = t / tu
    sxu = space_units_string(xu)
    syu = space_units_string(yu)
    szu = space_units_string(zu)
    stu = time_units_string(tu)

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    ix0 = argmin(abs.(x .- x0))
    iy0 = argmin(abs.(y .- y0))
    iz0 = argmin(abs.(z .- z0))

    isnothing(xlims[1]) ? xmin=x[1] : xmin=xlims[1]
    isnothing(xlims[2]) ? xmax=x[end] : xmax=xlims[2]
    isnothing(ylims[1]) ? ymin=y[1] : ymin=ylims[1]
    isnothing(ylims[2]) ? ymax=y[end] : ymax=ylims[2]
    isnothing(zlims[1]) ? zmin=z[1] : zmin=zlims[1]
    isnothing(zlims[2]) ? zmax=z[end] : zmax=zlims[2]

    fig = mak.Figure(resolution=(1600,600), fontsize=14)
    ax1 = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="y ($syu)")
    mak.xlims!(ax1, (xmin,xmax))
    mak.ylims!(ax1, (ymin,ymax))
    ax2 = mak.Axis(fig[1,2]; xlabel="x ($sxu)", ylabel="z ($szu)")
    mak.xlims!(ax2, (xmin,xmax))
    mak.ylims!(ax2, (zmin,zmax))
    ax3 = mak.Axis(fig[1,3]; xlabel="y ($syu)", ylabel="z ($szu)")
    mak.xlims!(ax3, (ymin,ymax))
    mak.ylims!(ax3, (zmin,zmax))
    mak.display(fig)

    colormap = cmap
    colorrange = (vmin,vmax)

    it = 1
    hm1 = mak.heatmap!(ax1, x, y, F[:,:,iz0,it]; colormap, colorrange)
    hm2 = mak.heatmap!(ax2, x, z, F[:,iy0,:,it]; colormap, colorrange)
    hm3 = mak.heatmap!(ax3, y, z, F[ix0,:,:,it]; colormap, colorrange)

    ax2.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
    mak.Colorbar(fig[2,3], hm1; label=svar, vertical=false, flipaxis=false)

    if movie
        ext = splitext(fname)[end]
        fname_movie = replace(fname, ext => ".mp4")
        mak.record(fig, fname_movie, 1:length(t); framerate=12) do it
            hm1[3] = F[:,:,iz0,it]
            hm2[3] = F[:,iy0,:,it]
            hm3[3] = F[ix0,:,:,it]
            ax2.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end
    else
        sg = mak.SliderGrid(fig[2,1:2], (label="Time", range=1:length(t), startvalue=1))
        mak.on(sg.sliders[1].value) do it
            hm1[3] = F[:,:,iz0,it]
            hm2[3] = F[:,iy0,:,it]
            hm3[3] = F[ix0,:,:,it]
            ax2.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end
    end
    return nothing
end


function inspect3D_poynting(
    fname; xu=1, yu=1, zu=1, tu=1, vmin=0, vmax=1, norm=false,
    cmap=mak.Reverse(:Hiroshige), aspect=(1,1,1),
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
    if norm
        F .= F ./ maximum(F)
    end

    @. x = x / xu
    @. y = y / yu
    @. z = z / zu
    @. t = t / tu
    sxu = space_units_string(xu)
    syu = space_units_string(yu)
    szu = space_units_string(zu)
    stu = time_units_string(tu)

    fig = mak.Figure(resolution=(950,992), fontsize=14)
    ax = mak.Axis3(
        fig[1,1];
        aspect, perspectiveness=0,
        xlabel="x ($sxu)", ylabel="y ($syu)", zlabel="z ($szu)",
    )
    mak.display(fig)

    cmap = mak.to_colormap(cmap)
    cmap[1] = mak.RGBAf(0,0,0,0)

    it = 1
    img = mak.volume!(
        ax, x, y, z, F[:,:,:,it];
        colormap=cmap, colorrange=(vmin,vmax),
        algorithm=:absorption, absorption=4f0,
        # algorithm=:iso, isovalue=0.1,
    )
    mak.Colorbar(fig[1,2], img; label="|S|")
    ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)

    sg = mak.SliderGrid(fig[2,1], (label="Time", range=1:length(t), startvalue=1))
    mak.on(sg.sliders[1].value) do it
        img[4] = F[:,:,:,it]
        ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
    end
    return nothing
end


function inspect3D_poynting_xsec(
    fname; xu=1, yu=1, zu=1, tu=1, vmin=0, vmax=1, norm=false,
    cmap=mak.Reverse(:Hiroshige),
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
    if norm
        F .= F ./ maximum(F)
    end

    @. x = x / xu
    @. y = y / yu
    @. z = z / zu
    @. t = t / tu
    sxu = space_units_string(xu)
    syu = space_units_string(yu)
    szu = space_units_string(zu)
    stu = time_units_string(tu)

    ix, iy, iz = (div(length(p),2) for p in (x,y,z))

    fig = mak.Figure(resolution=(1600,600), fontsize=14)
    ax1 = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="y ($syu)")
    ax2 = mak.Axis(fig[1,2]; xlabel="x ($sxu)", ylabel="z ($szu)")
    ax3 = mak.Axis(fig[1,3]; xlabel="y ($syu)", ylabel="z ($szu)")
    mak.display(fig)

    colormap = cmap
    colorrange = (vmin,vmax)

    it = 1
    hm1 = mak.heatmap!(ax1, x, y, F[:,:,iz,it]; colormap, colorrange)
    hm2 = mak.heatmap!(ax2, x, z, F[:,iy,:,it]; colormap, colorrange)
    hm3 = mak.heatmap!(ax3, y, z, F[ix,:,:,it]; colormap, colorrange)

    ax2.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
    mak.Colorbar(fig[2,3], hm1; label="|S|", vertical=false, flipaxis=false)

    sg = mak.SliderGrid(fig[2,1:2], (label="Time", range=1:length(t), startvalue=1))
    mak.on(sg.sliders[1].value) do it
        hm1[3] = F[:,:,iz,it]
        hm2[3] = F[:,iy,:,it]
        hm3[3] = F[ix,:,:,it]
        ax2.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
    end
    return nothing
end


function plot3D_poynting_averaged(
    fname; xu=1, yu=1, zu=1, vmin=0, vmax=1, norm=false,
    cmap=mak.Reverse(:Hiroshige), aspect=(1,1,1),
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
    F = dropdims(sum(F; dims=4); dims=4) ./ length(t)

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    @. x = x / xu
    @. y = y / yu
    @. z = z / zu
    sxu = space_units_string(xu)
    syu = space_units_string(yu)
    szu = space_units_string(zu)

    fig = mak.Figure(resolution=(950,992), fontsize=14)
    ax = mak.Axis3(
        fig[1,1];
        aspect, perspectiveness=0,
        xlabel="x ($sxu)", ylabel="y ($syu)", zlabel="z ($szu)",
    )
    mak.display(fig)

    cmap = mak.to_colormap(cmap)
    cmap[1] = mak.RGBAf(0,0,0,0)

    img = mak.volume!(
        ax, x, y, z, F;
        colormap=cmap, colorrange=(vmin,vmax),
        algorithm=:absorption, absorption=4f0,
    )
    mak.Colorbar(fig[1,2], img; label="Time averaged |S|")
    return nothing
end


function plot3D_poynting_averaged_xsec(
    fname; xu=1, yu=1, zu=1, vmin=0, vmax=1, norm=false,
    cmap=mak.Reverse(:Hiroshige), new_window=false,
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
    F = dropdims(sum(F; dims=4); dims=4) ./ length(t)

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    @. x = x / xu
    @. y = y / yu
    @. z = z / zu
    sxu = space_units_string(xu)
    syu = space_units_string(yu)
    szu = space_units_string(zu)

    ix, iy, iz = (div(length(p),2) for p in (x,y,z))

    fig = mak.Figure(resolution=(1600,600), fontsize=14)
    if new_window
        mak.display(mak.Screen(), fig)
    end
    ax1 = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="y ($syu)")
    ax2 = mak.Axis(fig[1,2]; xlabel="x ($sxu)", ylabel="z ($szu)")
    ax3 = mak.Axis(fig[1,3]; xlabel="y ($syu)", ylabel="z ($szu)")
    mak.display(fig)

    colormap = cmap
    colorrange = (vmin,vmax)

    hm1 = mak.heatmap!(ax1, x, y, F[:,:,iz]; colormap, colorrange)
    hm2 = mak.heatmap!(ax2, x, z, F[:,iy,:]; colormap, colorrange)
    hm3 = mak.heatmap!(ax3, y, z, F[ix,:,:]; colormap, colorrange)
    mak.Colorbar(fig[2,3], hm1; label="Time averaged |S|", vertical=false, flipaxis=false)
    return nothing
end
