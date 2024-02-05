# ******************************************************************************************
# Visualize rough edges, surfaces, and volumes
# ******************************************************************************************
function rough_plot(x, R; xu=nothing, Ru=nothing, new_window=false)
    isnothing(xu) ? xu = units(x) : nothing
    isnothing(Ru) ? Ru = units(R) : nothing
    sxu = units_name_space(xu)
    sRu = units_name_space(Ru)

    fig = mak.Figure(size=(800,600))
    ax = mak.Axis(fig[1,1]; xlabel="x ($sxu)", ylabel="R ($sRu)")
    mak.lines!(ax, x/xu, R/Ru)
    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end
    return nothing
end


function rough_plot(
    x, y, R;
    xu=nothing, yu=nothing, Ru=nothing, colormap=CMAPDIV, colorrange=nothing, save=false,
    save_fname="out.png", new_window=false,
)
    isnothing(xu) ? xu = units(x) : nothing
    isnothing(yu) ? yu = units(y) : nothing
    isnothing(Ru) ? Ru = units(R) : nothing
    sxu = units_name_space(xu)
    syu = units_name_space(yu)
    sRu = units_name_space(Ru)

    if isnothing(colorrange)
        Rmax = maximum(abs, R/Ru)
        colorrange = (-Rmax, Rmax)
    end

    x1, x2 = x[1]/xu, x[end]/xu
    y1, y2 = y[1]/yu, y[end]/yu

    Nx, Ny = size(R)
    Rx = R[:,halfint(Ny)]
    Ry = R[halfint(Nx),:]

    fig = mak.Figure(size=(950,992))
    ax11 = mak.Axis(fig[1,1]; ylabel="R ($sRu)")
    ax21 = mak.Axis(fig[2,1]; xlabel="x ($sxu)", ylabel="y ($syu)")
    ax22 = mak.Axis(fig[2,2]; xlabel="R ($sRu)")

    mak.xlims!(ax11, (x1,x2))
    mak.ylims!(ax11, colorrange...)
    mak.xlims!(ax22, colorrange...)
    mak.ylims!(ax22, (y1,y2))

    mak.linkxaxes!(ax11, ax21)
    mak.linkyaxes!(ax22, ax21)
    mak.hidexdecorations!(ax11, ticks=false)
    mak.hideydecorations!(ax22, ticks=false)

    mak.rowsize!(fig.layout, 1, mak.Relative(1/4))
    mak.colsize!(fig.layout, 1, mak.Relative(3/4))

    mak.lines!(ax11, x/xu, Rx/Ru)

    hm = mak.heatmap!(ax21, x/xu, y/yu, R/Ru; colormap, colorrange)
    mak.lines!(ax21, [x1,x2], [0,0]; color=:black, linewidth=1)
    mak.lines!(ax21, [0,0], [y1,y2]; color=:black, linewidth=1)
    mak.Colorbar(fig[3,1], hm; vertical=false, flipaxis=false, label="R ($sRu)")

    mak.lines!(ax22, Ry/Ru, y/yu)

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


function rough_plot(
    x, y, z, R;
    xu=nothing, yu=nothing, zu=nothing, Ru=nothing, colormap=CMAPDIV, colorrange=nothing,
    aspect=(1,1,1), colorbar=true, algorithm=:volumeslices, absorption=1, isovalue=0,
    new_window=false,
)
    isnothing(xu) ? xu = units(x) : nothing
    isnothing(yu) ? yu = units(y) : nothing
    isnothing(zu) ? zu = units(z) : nothing
    isnothing(Ru) ? Ru = units(R) : nothing
    sxu = units_name_space(xu)
    syu = units_name_space(yu)
    szu = units_name_space(zu)
    sRu = units_name_space(Ru)

    if isnothing(colorrange)
        Rmax = maximum(abs, R/Ru)
        colorrange = (-Rmax, Rmax)
    end

    Nx, Ny, Nz = size(R)
    ix, iy, iz = halfint(Nx), halfint(Ny), halfint(Nz)

    fig = mak.Figure(size=(950,992))
    ax = mak.Axis3(
        fig[1,1];
        xlabel="x ($sxu)", ylabel="y ($syu)", zlabel="z ($szu)", aspect, perspectiveness=0,
    )

    if algorithm == :volumeslices
        vol = mak.volumeslices!(
            ax, x/xu, y/yu, z/zu, R/Ru; bbox_visible=false, colormap, colorrange,
        )
        vol.update_yz[](ix)
        vol.update_xz[](iy)
        vol.update_xy[](iz)

        if colorbar
            mak.Colorbar(
                fig[1,2], vol.heatmap_xy[], height=mak.Relative(0.5), label="R ($sRu)",
            )
        end

        sg = mak.SliderGrid(
            fig[2,1],
            (label="x", range=1:length(x), startvalue=ix),
            (label="y", range=1:length(y), startvalue=iy),
            (label="z", range=1:length(z), startvalue=iz),
        )
        mak.on(sg.sliders[1].value) do i
            vol.update_yz[](i)
        end
        mak.on(sg.sliders[2].value) do i
            vol.update_xz[](i)
        end
        mak.on(sg.sliders[3].value) do i
            vol.update_xy[](i)
        end

        # hmaps = [vol.heatmap_yz[], vol.heatmap_xz[], vol.heatmap_xy[]]
        # toggles = [mak.Toggle(sg.layout[i,4], active=true) for i in 1:length(hmaps)]
        # map(zip(hmaps, toggles)) do (h, t)
        #     mak.connect!(h.visible, t.active)
        # end
    else
        Rmax = maximum(abs, R/Ru)

        img = mak.volume!(
            ax, x/xu, y/yu, z/zu, R/Ru;
            colormap, colorrange, algorithm, absorption, isovalue, isorange=0.05*Rmax,
        )

        if colorbar
            mak.Colorbar(fig[1,2], img, height=mak.Relative(0.5))
        end
    end

    if new_window
        mak.display(mak.Screen(), fig)
    else
        mak.display(fig)
    end

    return nothing
end
