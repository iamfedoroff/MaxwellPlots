function inspect(
    x, y, z, t, F;
    xu=1, yu=1, zu=1, tu=1, norm=true, colormap=nothing, colorrange=(nothing,nothing),
    aspect=:data, xlims=(nothing,nothing), ylims=(nothing,nothing), zlims=(nothing,nothing),
    tlims=(nothing,nothing), new_window=false, movie=false, movie_fname="movie.mp4",
    colorbar=true,
)
    @. x = x / xu
    @. y = y / yu
    @. z = z / zu
    @. t = t / tu
    sxu = space_units_string(xu)
    syu = space_units_string(yu)
    szu = space_units_string(zu)
    stu = time_units_string(tu)

    x, y, z, t, F = apply_limits(x, y, z, t, F; xlims, ylims, zlims, tlims)

    it = 1

    @show extrema(F)
    if norm
        F .= F ./ maximum(F)
    end


    # --------------------------------------------------------------------------------------
    vmin, vmax = colorrange
    isnothing(vmin) ? vmin = minimum(F) : nothing
    isnothing(vmax) ? vmax = maximum(F) : nothing
    colorrange = (vmin, vmax)
    if vmin * vmax < 0
        isnothing(colormap) ? colormap = :seismic : nothing
        colormap = mak.to_colormap(colormap)
        Nmid = halfint(length(colormap))
        @. colormap[Nmid-1:Nmid+1] = mak.RGBAf(0,0,0,0)
    else
        isnothing(colormap) ? colormap = mak.Reverse(:Hiroshige) : nothing
        colormap = mak.to_colormap(colormap)
        colormap[1] = mak.RGBAf(0,0,0,0)
    end

    fig = mak.Figure(resolution=(950,992), fontsize=14)

    ax = mak.Axis3(
        fig[1,1];
        aspect, perspectiveness=0, xlabel="x ($sxu)", ylabel="y ($syu)", zlabel="z ($szu)",
    )
    mak.limits!(ax, x[1], x[end], y[1], y[end], z[1], z[end])
    ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)

    img = mak.volume!(
        ax, x, y, z, F[:,:,:,it]; colormap, colorrange,
        algorithm=:absorption, absorption=4f0,
        # algorithm=:iso, isovalue=0.5*vmax,
    )
    # img = mak.contour!(
    #     ax, x, y, z, F[:,:,:,it];
    #     levels=[0.1*vmin,0.1*vmax], colormap=cmap, colorrange=(vmin,vmax),
    #     alpha=1,
    # )

    if colorbar
        mak.Colorbar(fig[1,2], img; height=mak.Relative(0.5))
    end

    if movie
        mak.record(fig, movie_fname, 1:length(t); framerate=6) do it
            img[4] = F[:,:,:,it]
            ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end
    else
        sg = mak.SliderGrid(fig[2,1], (label="Time", range=1:length(t), startvalue=1))
        mak.on(sg.sliders[1].value) do it
            img[4] = F[:,:,:,it]
            ax.title[] = @sprintf("%d:     %.3f (%s)", it, t[it], stu)
        end
    end

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end
    return nothing
end


function inspect_xsec(
    x, y, z, t, F, x0, y0, z0;
    xu, yu, zu, tu, norm, norm_point, vmin, vmax, aspect, xlims, ylims, zlims, cmap,
    new_window, movie, movie_fname, guidelines=true,
)
    @. x = x / xu
    @. y = y / yu
    @. z = z / zu
    @. t = t / tu
    sxu = space_units_string(xu)
    syu = space_units_string(yu)
    szu = space_units_string(zu)
    stu = time_units_string(tu)

    ix0 = argmin(abs.(x .- x0))
    iy0 = argmin(abs.(y .- y0))
    iz0 = argmin(abs.(z .- z0))

    if norm
        if isnothing(norm_point)
            F .= F ./ maximum(F)
        else
            xn, yn, zn = norm_point
            ixn = argmin(abs.(x .- xn))
            iyn = argmin(abs.(y .- yn))
            izn = argmin(abs.(z .- zn))
            @views F .= F ./ maximum(F[ixn,iyn,izn,:])
        end
    end

    isnothing(xlims[1]) ? xmin=x[1] : xmin=xlims[1]
    isnothing(xlims[2]) ? xmax=x[end] : xmax=xlims[2]
    isnothing(ylims[1]) ? ymin=y[1] : ymin=ylims[1]
    isnothing(ylims[2]) ? ymax=y[end] : ymax=ylims[2]
    isnothing(zlims[1]) ? zmin=z[1] : zmin=zlims[1]
    isnothing(zlims[2]) ? zmax=z[end] : zmax=zlims[2]

    fig = mak.Figure(resolution=(1600,600), fontsize=14)
    ax1 = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="y ($syu)", aspect=aspect[1])
    ax2 = mak.Axis(fig[1,2]; xlabel="x ($sxu)", ylabel="z ($szu)", aspect=aspect[2])
    ax3 = mak.Axis(fig[1,3]; xlabel="y ($syu)", ylabel="z ($szu)", aspect=aspect[3])
    mak.xlims!(ax1, (xmin,xmax))
    mak.ylims!(ax1, (ymin,ymax))
    mak.xlims!(ax2, (xmin,xmax))
    mak.ylims!(ax2, (zmin,zmax))
    mak.xlims!(ax3, (ymin,ymax))
    mak.ylims!(ax3, (zmin,zmax))
    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end

    colormap = cmap
    colorrange = (vmin,vmax)

    it = 1
    hm1 = mak.heatmap!(ax1, x, y, F[:,:,iz0,it]; colormap, colorrange)
    hm2 = mak.heatmap!(ax2, x, z, F[:,iy0,:,it]; colormap, colorrange)
    hm3 = mak.heatmap!(ax3, y, z, F[ix0,:,:,it]; colormap, colorrange)
    mak.Colorbar(fig[2,3], hm1; vertical=false, flipaxis=false)

    # guide lines:
    if guidelines
        x0 = x[ix0]
        y0 = y[iy0]
        z0 = z[iz0]
        mak.lines!(ax1, [xmin, xmax], [y0, y0]; color=:white, linewidth=0.5)
        mak.lines!(ax1, [x0, x0], [ymin, ymax]; color=:white, linewidth=0.5)
        mak.lines!(ax2, [xmin, xmax], [z0, z0]; color=:white, linewidth=0.5)
        mak.lines!(ax2, [x0, x0], [zmin, zmax]; color=:white, linewidth=0.5)
        mak.lines!(ax3, [ymin, ymax], [z0, z0]; color=:white, linewidth=0.5)
        mak.lines!(ax3, [y0, y0], [zmin, zmax]; color=:white, linewidth=0.5)
    end

    if movie
        mak.record(fig, movie_fname, 1:length(t); framerate=12) do it
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


function mplot(
    x, y, z, F;
    xu, yu, zu, norm, vmin, vmax, aspect, xlims, ylims, zlims, cmap, new_window, save,
)
    @. x = x / xu
    @. y = y / yu
    @. z = z / zu
    sxu = space_units_string(xu)
    syu = space_units_string(yu)
    szu = space_units_string(zu)

    if norm
        F .= F ./ maximum(F)
    end

    isnothing(xlims[1]) ? xmin=x[1] : xmin=xlims[1]
    isnothing(xlims[2]) ? xmax=x[end] : xmax=xlims[2]
    isnothing(ylims[1]) ? ymin=y[1] : ymin=ylims[1]
    isnothing(ylims[2]) ? ymax=y[end] : ymax=ylims[2]
    isnothing(zlims[1]) ? zmin=z[1] : zmin=zlims[1]
    isnothing(zlims[2]) ? zmax=z[end] : zmax=zlims[2]

    fig = mak.Figure(resolution=(950,992), fontsize=14)
    ax = mak.Axis3(
        fig[1,1];
        aspect, perspectiveness=0, xlabel="x ($sxu)", ylabel="y ($syu)", zlabel="z ($szu)",
    )
    mak.xlims!(ax, (xmin,xmax))
    mak.ylims!(ax, (ymin,ymax))
    mak.zlims!(ax, (zmin,zmax))
    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end

    if vmin * vmax < 0
        colormap = mak.to_colormap(cmap)
        Nmid = div(length(colormap),2)
        # colormap[Nmid+1] = mak.RGBAf(0,0,0,0)
        @. colormap[Nmid:Nmid+2] = mak.RGBAf(0,0,0,0)
    else
        colormap = mak.to_colormap(cmap)
        colormap[1] = mak.RGBAf(0,0,0,0)
    end

    img = mak.volume!(
        ax, x, y, z, F;
        colormap, colorrange=(vmin,vmax),
        algorithm=:absorption, absorption=4f0,
    )
    mak.Colorbar(fig[1,2], img)

    if save
        ext = splitext(fname)[end]
        fname_fig = replace(fname, ext => ".png")
        mak.save(fname_fig, fig)
    end
    return nothing
end


function mplot_xsec(
    x, y, z, F, x0, y0, z0;
    xu, yu, zu, norm, norm_point, vmin, vmax, aspect, xlims, ylims, zlims, cmap,
    new_window, save, fname_fig, guidelines=true,
)

    @. x = x / xu
    @. y = y / yu
    @. z = z / zu
    sxu = space_units_string(xu)
    syu = space_units_string(yu)
    szu = space_units_string(zu)

    ix0 = argmin(abs.(x .- x0))
    iy0 = argmin(abs.(y .- y0))
    iz0 = argmin(abs.(z .- z0))

    if norm
        if isnothing(norm_point)
            F .= F ./ maximum(F)
        else
            xn, yn, zn = norm_point
            ixn = argmin(abs.(x .- xn))
            iyn = argmin(abs.(y .- yn))
            izn = argmin(abs.(z .- zn))
            @views F .= F ./ maximum(F[ixn,iyn,izn,:])
        end
    end

    isnothing(xlims[1]) ? xmin=x[1] : xmin=xlims[1]
    isnothing(xlims[2]) ? xmax=x[end] : xmax=xlims[2]
    isnothing(ylims[1]) ? ymin=y[1] : ymin=ylims[1]
    isnothing(ylims[2]) ? ymax=y[end] : ymax=ylims[2]
    isnothing(zlims[1]) ? zmin=z[1] : zmin=zlims[1]
    isnothing(zlims[2]) ? zmax=z[end] : zmax=zlims[2]

    fig = mak.Figure(resolution=(1600,600), fontsize=14)
    ax1 = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="y ($syu)", aspect=aspect[1])
    ax2 = mak.Axis(fig[1,2]; xlabel="x ($sxu)", ylabel="z ($szu)", aspect=aspect[2])
    ax3 = mak.Axis(fig[1,3]; xlabel="y ($syu)", ylabel="z ($szu)", aspect=aspect[3])
    mak.xlims!(ax1, (xmin,xmax))
    mak.ylims!(ax1, (ymin,ymax))
    mak.xlims!(ax2, (xmin,xmax))
    mak.ylims!(ax2, (zmin,zmax))
    mak.xlims!(ax3, (ymin,ymax))
    mak.ylims!(ax3, (zmin,zmax))
    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end

    colormap = cmap
    colorrange = (vmin,vmax)

    hm1 = mak.heatmap!(ax1, x, y, F[:,:,iz0]; colormap, colorrange)
    hm2 = mak.heatmap!(ax2, x, z, F[:,iy0,:]; colormap, colorrange)
    hm3 = mak.heatmap!(ax3, y, z, F[ix0,:,:]; colormap, colorrange)
    mak.Colorbar(fig[2,3], hm1; vertical=false, flipaxis=false)

    # guide lines:
    if guidelines
        x0 = x[ix0]
        y0 = y[iy0]
        z0 = z[iz0]
        mak.lines!(ax1, [xmin, xmax], [y0, y0]; color=:white, linewidth=0.5)
        mak.lines!(ax1, [x0, x0], [ymin, ymax]; color=:white, linewidth=0.5)
        mak.lines!(ax2, [xmin, xmax], [z0, z0]; color=:white, linewidth=0.5)
        mak.lines!(ax2, [x0, x0], [zmin, zmax]; color=:white, linewidth=0.5)
        mak.lines!(ax3, [ymin, ymax], [z0, z0]; color=:white, linewidth=0.5)
        mak.lines!(ax3, [y0, y0], [zmin, zmax]; color=:white, linewidth=0.5)
    end

    if save
        mak.save(fname_fig, fig)
    end
    return nothing
end
