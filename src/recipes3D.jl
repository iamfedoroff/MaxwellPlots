function inspect_xsec(
    x, y, z, t, F;
    xu=1, yu=1, zu=1, tu=1,
    xcut=nothing, ycut=nothing, zcut=nothing,
    xlims=nothing, ylims=nothing, zlims=nothing, tlims=nothing,
    norm=true, colormap=nothing, colorrange=nothing, colorbar=true, aspect=(1,1,1),
    movie=true, movie_fname="movie.mp4", new_window=false, guidelines=true,
)
    x, y, z, t, F = apply_limits(x, y, z, t, F; xlims, ylims, zlims, tlims)

    Nx, Ny, Nz, Nt = size(F)
    isnothing(xcut) ? ix = halfint(Nx) : ix = argmin(abs.(x.-xcut))
    isnothing(ycut) ? iy = halfint(Ny) : iy = argmin(abs.(y.-ycut))
    isnothing(zcut) ? iz = halfint(Nz) : iz = argmin(abs.(z.-zcut))
    Lx, Ly, Lz  = x[end]-x[1], y[end]-y[1], z[end]-z[1]

    xu = isnothing(xu) ? units(x) : xu
    yu = isnothing(yu) ? units(y) : yu
    zu = isnothing(zu) ? units(z) : zu
    tu = isnothing(tu) ? units(t) : tu
    sxu = units_name_space(xu)
    syu = units_name_space(yu)
    szu = units_name_space(zu)
    stu = units_name_time(tu)

    @. x = x / xu
    @. y = y / yu
    @. z = z / zu
    @. t = t / tu

    it = 1

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    # --------------------------------------------------------------------------------------
    isnothing(colorrange) ? colorrange = (minimum(F), maximum(F)) : nothing
    vmin, vmax = colorrange
    isnothing(vmin) ? vmin = minimum(F) : nothing
    isnothing(vmax) ? vmax = maximum(F) : nothing
    if vmin * vmax < 0   # diverging colormap
        isnothing(colormap) ? colormap = CMAPDIV : nothing
    else
        isnothing(colormap) ? colormap = CMAP : nothing
    end

    isnothing(aspect) ? aspect = (Lx/Ly, Lx/Lz, Ly/Lz) : nothing
    isnothing(aspect[1]) ? aspect[1] = Lx/Ly : nothing
    isnothing(aspect[2]) ? aspect[2] = Lx/Lz : nothing
    isnothing(aspect[3]) ? aspect[3] = Ly/Lz : nothing

    fig = mak.Figure(size=(1600,600))
    ax1 = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="y ($syu)", aspect=aspect[1])
    ax2 = mak.Axis(fig[1,2]; xlabel="x ($sxu)", ylabel="z ($szu)", aspect=aspect[2])
    ax3 = mak.Axis(fig[1,3]; xlabel="y ($syu)", ylabel="z ($szu)", aspect=aspect[3])
    mak.limits!(ax1, x[1], x[end], y[1], y[end])
    mak.limits!(ax2, x[1], x[end], z[1], z[end])
    mak.limits!(ax3, y[1], y[end], z[1], z[end])
    ax2.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)

    hm1 = mak.heatmap!(ax1, x, y, F[:,:,iz,it]; colormap, colorrange)
    hm2 = mak.heatmap!(ax2, x, z, F[:,iy,:,it]; colormap, colorrange)
    hm3 = mak.heatmap!(ax3, y, z, F[ix,:,:,it]; colormap, colorrange)

    if colorbar
        mak.Colorbar(fig[2,3], hm1; vertical=false, flipaxis=false)
    end

    # guide lines:
    if guidelines
        color = :gray15
        linewidth = 0.5
        mak.lines!(ax1, [x[1], x[end]], [y[iy], y[iy]]; color, linewidth)
        mak.lines!(ax1, [x[ix], x[ix]], [y[1], y[end]]; color, linewidth)
        mak.lines!(ax2, [x[1], x[end]], [z[iz], z[iz]]; color, linewidth)
        mak.lines!(ax2, [x[ix], x[ix]], [z[1], z[end]]; color, linewidth)
        mak.lines!(ax3, [y[1], y[end]], [z[iz], z[iz]]; color, linewidth)
        mak.lines!(ax3, [y[iy], y[iy]], [z[1], z[end]]; color, linewidth)
    end

    if movie
        mak.record(fig, movie_fname, 1:Nt; framerate=12) do it
            hm1[3] = F[:,:,iz,it]
            hm2[3] = F[:,iy,:,it]
            hm3[3] = F[ix,:,:,it]
            ax2.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end
    else
        sg = mak.SliderGrid(fig[2,1:2], (label="Time", range=1:Nt, startvalue=1))
        mak.on(sg.sliders[1].value) do it
            hm1[3] = F[:,:,iz,it]
            hm2[3] = F[:,iy,:,it]
            hm3[3] = F[ix,:,:,it]
            ax2.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end
    end

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end
    return nothing
end


function inspect_volume(
    x, y, z, F;
    xu=nothing, yu=nothing, zu=nothing,
    xcut=nothing, ycut=nothing, zcut=nothing,
    xlims=nothing, ylims=nothing, zlims=nothing,
    norm=true, colormap=nothing, colorrange=nothing, colorbar=true, aspect=(1,1,1),
    save=false, save_fname="image.png", new_window=false,
)
    x, y, z, F = apply_limits(x, y, z, F; xlims, ylims, zlims)

    Nx, Ny, Nz = size(F)
    isnothing(xcut) ? ix = halfint(Nx) : ix = argmin(abs.(x.-xcut))
    isnothing(ycut) ? iy = halfint(Ny) : iy = argmin(abs.(y.-ycut))
    isnothing(zcut) ? iz = halfint(Nz) : iz = argmin(abs.(z.-zcut))

    xu = isnothing(xu) ? units(x) : xu
    yu = isnothing(yu) ? units(y) : yu
    zu = isnothing(zu) ? units(z) : zu
    sxu = units_name_space(xu)
    syu = units_name_space(yu)
    szu = units_name_space(zu)

    @. x = x / xu
    @. y = y / yu
    @. z = z / zu

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    # --------------------------------------------------------------------------------------
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

    ax = mak.Axis3(
        fig[1,1]; xlabel="x ($sxu)", ylabel="y ($syu)", zlabel="z ($szu)", aspect,
        perspectiveness=0,
    )
    mak.limits!(ax, x[1], x[end], y[1], y[end], z[1], z[end])
    ax.title[] = @sprintf(
        "x=%.3f %s   y=%.3f %s   z=%.3f %s", x[ix], sxu, y[iy], syu, z[iz], szu,
    )

    vol = mak.volumeslices!(
        ax, x, y, z, F; bbox_visible=false, colormap, colorrange,
    )
    vol.update_yz[](ix)
    vol.update_xz[](iy)
    vol.update_xy[](iz)

    # Spines for cut planes:
    # color = :gray15
    # linewidth = 0.7
    # ix > 1  ? ixm1=ix-1 : ixm1=1
    # ix < Nx ? ixp1=ix+1 : ixp1=Nx
    # iy > 1  ? iym1=iy-1 : iym1=1
    # iy < Ny ? iyp1=iy+1 : iyp1=Ny
    # iz > 1  ? izm1=iz-1 : izm1=1
    # iz < Nz ? izp1=iz+1 : izp1=Nz
    # # x
    # mak.lines!(ax, [x[1],x[end]], [y[iym1],y[iym1]], [z[izm1],z[izm1]]; color, linewidth)
    # mak.lines!(ax, [x[1],x[end]], [y[iym1],y[iym1]], [z[izp1],z[izp1]]; color, linewidth)
    # mak.lines!(ax, [x[1],x[end]], [y[iyp1],y[iyp1]], [z[izm1],z[izm1]]; color, linewidth)
    # mak.lines!(ax, [x[1],x[end]], [y[iyp1],y[iyp1]], [z[izp1],z[izp1]]; color, linewidth)
    # # y
    # mak.lines!(ax, [x[ixm1],x[ixm1]], [y[1],y[end]], [z[izm1],z[izm1]]; color, linewidth)
    # mak.lines!(ax, [x[ixm1],x[ixm1]], [y[1],y[end]], [z[izp1],z[izp1]]; color, linewidth)
    # mak.lines!(ax, [x[ixp1],x[ixp1]], [y[1],y[end]], [z[izm1],z[izm1]]; color, linewidth)
    # mak.lines!(ax, [x[ixp1],x[ixp1]], [y[1],y[end]], [z[izp1],z[izp1]]; color, linewidth)
    # # z
    # mak.lines!(ax, [x[ixm1],x[ixm1]], [y[iym1],y[iym1]], [z[1],z[end]]; color, linewidth)
    # mak.lines!(ax, [x[ixm1],x[ixm1]], [y[iyp1],y[iyp1]], [z[1],z[end]]; color, linewidth)
    # mak.lines!(ax, [x[ixp1],x[ixp1]], [y[iym1],y[iym1]], [z[1],z[end]]; color, linewidth)
    # mak.lines!(ax, [x[ixp1],x[ixp1]], [y[iyp1],y[iyp1]], [z[1],z[end]]; color, linewidth)
    # mak.hidedecorations!(ax)
    # mak.hidespines!(ax)

    if colorbar
        mak.Colorbar(fig[1,2], vol.heatmap_xy[], height=mak.Relative(0.5))
    end

    if save
        mak.save(save_fname, fig)
    else
        sg = mak.SliderGrid(
            fig[2,1],
            (label="x", range=1:length(x), startvalue=ix),
            (label="y", range=1:length(y), startvalue=iy),
            (label="z", range=1:length(z), startvalue=iz),
        )
        mak.on(sg.sliders[1].value) do i
            ix = i
            vol.update_yz[](ix)
            ax.title[] = @sprintf(
                "x=%.3f %s   y=%.3f %s   z=%.3f %s", x[ix], sxu, y[iy], syu, z[iz], szu,
            )
        end
        mak.on(sg.sliders[2].value) do i
            iy = i
            vol.update_xz[](iy)
            ax.title[] = @sprintf(
                "x=%.3f %s   y=%.3f %s   z=%.3f %s", x[ix], sxu, y[iy], syu, z[iz], szu,
            )
        end
        mak.on(sg.sliders[3].value) do i
            iz = i
            vol.update_xy[](iz)
            ax.title[] = @sprintf(
                "x=%.3f %s   y=%.3f %s   z=%.3f %s", x[ix], sxu, y[iy], syu, z[iz], szu,
            )
        end

        # hmaps = [vol.heatmap_yz[], vol.heatmap_xz[], vol.heatmap_xy[]]
        # toggles = [mak.Toggle(sg.layout[i,4], active=true) for i in 1:length(hmaps)]
        # map(zip(hmaps, toggles)) do (h, t)
        #     mak.connect!(h.visible, t.active)
        # end
    end

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end
    return nothing
end


function plot_volume(
    x, y, z, F;
    xu=nothing, yu=nothing, zu=nothing,
    xlims=nothing, ylims=nothing, zlims=nothing,
    norm=true, colormap=nothing, colorrange=nothing, colorbar=true, aspect=(1,1,1),
    save=false, save_fname="image.png", new_window=false,
)
    x, y, z, F = apply_limits(x, y, z, F; xlims, ylims, zlims)

    xu = isnothing(xu) ? units(x) : xu
    yu = isnothing(yu) ? units(y) : yu
    zu = isnothing(zu) ? units(z) : zu
    sxu = units_name_space(xu)
    syu = units_name_space(yu)
    szu = units_name_space(zu)

    @. x = x / xu
    @. y = y / yu
    @. z = z / zu

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    # --------------------------------------------------------------------------------------
    isnothing(colorrange) ? colorrange = (minimum(F), maximum(F)) : nothing
    vmin, vmax = colorrange
    isnothing(vmin) ? vmin = minimum(F) : nothing
    isnothing(vmax) ? vmax = maximum(F) : nothing
    colorrange = (vmin, vmax)
    if vmin * vmax < 0
        isnothing(colormap) ? colormap = CMAPDIV : nothing
        colormap = mak.to_colormap(colormap)
        Nmid = halfint(length(colormap))
        @. colormap[Nmid-1:Nmid+1] = mak.RGBAf(0,0,0,0)
    else
        isnothing(colormap) ? colormap = CMAP : nothing
        colormap = mak.to_colormap(colormap)
        colormap[1] = mak.RGBAf(0,0,0,0)
    end

    fig = mak.Figure(size=(950,992))

    ax = mak.Axis3(
        fig[1,1]; xlabel="x ($sxu)", ylabel="y ($syu)", zlabel="z ($szu)", aspect,
        perspectiveness=0,
    )
    mak.limits!(ax, x[1], x[end], y[1], y[end], z[1], z[end])

    xrange = (x[1], x[end])
    yrange = (y[1], y[end])
    zrange = (z[1], z[end])

    img = mak.volume!(
        ax, xrange, yrange, zrange, F; colormap, colorrange,
        algorithm=:absorption, absorption=4f0,
    )

    if colorbar
        mak.Colorbar(fig[1,2], img; height=mak.Relative(0.5))
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


function plot_volume_xsec(
    x, y, z, F;
    xu=nothing, yu=nothing, zu=nothing,
    xcut=nothing, ycut=nothing, zcut=nothing,
    xlims=nothing, ylims=nothing, zlims=nothing,
    norm=true, colormap=nothing, colorrange=nothing, colorbar=true, aspect=(1,1,1),
    save=false, save_fname="image.png", new_window=false, guidelines=true,
)
    x, y, z, F = apply_limits(x, y, z, F; xlims, ylims, zlims)

    Nx, Ny, Nz = size(F)
    isnothing(xcut) ? ix = halfint(Nx) : ix = argmin(abs.(x.-xcut))
    isnothing(ycut) ? iy = halfint(Ny) : iy = argmin(abs.(y.-ycut))
    isnothing(zcut) ? iz = halfint(Nz) : iz = argmin(abs.(z.-zcut))
    Lx, Ly, Lz  = x[end]-x[1], y[end]-y[1], z[end]-z[1]

    xu = isnothing(xu) ? units(x) : xu
    yu = isnothing(yu) ? units(y) : yu
    zu = isnothing(zu) ? units(z) : zu
    sxu = units_name_space(xu)
    syu = units_name_space(yu)
    szu = units_name_space(zu)

    @. x = x / xu
    @. y = y / yu
    @. z = z / zu

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end

    # --------------------------------------------------------------------------------------
    isnothing(colorrange) ? colorrange = (minimum(F), maximum(F)) : nothing
    vmin, vmax = colorrange
    isnothing(vmin) ? vmin = minimum(F) : nothing
    isnothing(vmax) ? vmax = maximum(F) : nothing
    if vmin * vmax < 0   # diverging colormap
        isnothing(colormap) ? colormap = CMAPDIV : nothing
    else
        isnothing(colormap) ? colormap = CMAP : nothing
    end

    isnothing(aspect) ? aspect = (Lx/Ly, Lx/Lz, Ly/Lz) : nothing
    isnothing(aspect[1]) ? aspect[1] = Lx/Ly : nothing
    isnothing(aspect[2]) ? aspect[2] = Lx/Lz : nothing
    isnothing(aspect[3]) ? aspect[3] = Ly/Lz : nothing

    fig = mak.Figure(size=(1600,600))
    ax1 = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="y ($syu)", aspect=aspect[1])
    ax2 = mak.Axis(fig[1,2]; xlabel="x ($sxu)", ylabel="z ($szu)", aspect=aspect[2])
    ax3 = mak.Axis(fig[1,3]; xlabel="y ($syu)", ylabel="z ($szu)", aspect=aspect[3])
    mak.limits!(ax1, x[1], x[end], y[1], y[end])
    mak.limits!(ax2, x[1], x[end], z[1], z[end])
    mak.limits!(ax3, y[1], y[end], z[1], z[end])

    hm1 = mak.heatmap!(ax1, x, y, F[:,:,iz]; colormap, colorrange)
    hm2 = mak.heatmap!(ax2, x, z, F[:,iy,:]; colormap, colorrange)
    hm3 = mak.heatmap!(ax3, y, z, F[ix,:,:]; colormap, colorrange)

    if colorbar
        mak.Colorbar(fig[2,2], hm1; vertical=false, flipaxis=false)
    end

    # guide lines:
    if guidelines
        color = :gray15
        linewidth = 0.5
        mak.lines!(ax1, [x[1], x[end]], [y[iy], y[iy]]; color, linewidth)
        mak.lines!(ax1, [x[ix], x[ix]], [y[1], y[end]]; color, linewidth)
        mak.lines!(ax2, [x[1], x[end]], [z[iz], z[iz]]; color, linewidth)
        mak.lines!(ax2, [x[ix], x[ix]], [z[1], z[end]]; color, linewidth)
        mak.lines!(ax3, [y[1], y[end]], [z[iz], z[iz]]; color, linewidth)
        mak.lines!(ax3, [y[iy], y[iy]], [z[1], z[end]]; color, linewidth)
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
