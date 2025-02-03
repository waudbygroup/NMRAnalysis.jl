function gui!(state, expt::FixedPeakExperiment)
    g = Dict{Symbol, Any}() # GUI state
    state[:gui] = g

    g[:fig] = Figure(size=(1000,700))
    # g[:paneltop] = g[:fig][1,1] = GridLayout()
    # g[:panelmiddle] = g[:fig][2,1] = GridLayout()
    # g[:panelcontour] = g[:panelmiddle][1,1] = GridLayout()
    # # g[:panel3d] = g[:panelmiddle][1,2] = GridLayout()
    # g[:panelbottom] = g[:fig][3,1] = GridLayout()
    # g[:panelinfo] = g[:panelbottom][1,1] = GridLayout()
    # g[:panelpeakplot] = g[:panelbottom][1,2] = GridLayout()
    # rowsize!(g[:fig].layout, 1, Auto(true))
    # rowsize!(g[:fig].layout, 2, Auto(2))
    # rowsize!(g[:fig].layout, 3, Auto(1))
    g[:paneltop] = g[:fig][1,1:2] = GridLayout()
    g[:panelcontour] = g[:fig][2:3,1] = GridLayout()
    g[:panelinfo] = g[:fig][2,2] = GridLayout()
    g[:panelpeakplot] = g[:fig][3,2] = GridLayout()
    rowsize!(g[:fig].layout, 1, Auto(true))
    rowsize!(g[:fig].layout, 2, Auto(false,1))
    rowsize!(g[:fig].layout, 3, Auto(false,1))
    colsize!(g[:fig].layout, 1, Auto(false,2))
    colsize!(g[:fig].layout, 2, Auto(false,1))


    # top panel
    g[:cmdcontourup] = Button(g[:paneltop][1,1], label="contour ↑")
    g[:cmdcontourdown] = Button(g[:paneltop][1,2], label="contour ↓")
    g[:sliderslice] = Slider(g[:paneltop][1,3], range = 1:nslices(expt))
    g[:cmdsliceleft] = Button(g[:paneltop][1,4], label="←")
    g[:cmdsliceright] = Button(g[:paneltop][1,5], label="→")
    g[:slicelabel] = Label(g[:paneltop][1,6], state[:current_slice_label])

    # create contour plot
    g[:basecontour] = Observable(10.0)
    g[:contourscale] = Observable(1.7)
    g[:logscale] = lift(r -> r .^ (0:10), g[:contourscale])
    g[:contourlevels] = lift((c0,arr) -> c0 * arr, g[:basecontour], g[:logscale])

    g[:axcontour] = Axis(g[:panelcontour][1,1], xlabel="δX / ppm", ylabel="δy / ppm",
        xreversed=true, yreversed=true)
    g[:pltmask] = heatmap!(g[:axcontour],
        state[:current_mask_x],
        state[:current_mask_y],
        state[:current_mask_z],
        colormap=[:white,:lightgoldenrod1],
        colorrange=(0,1),
        )
    g[:pltcontour] = contour!(g[:axcontour],
        state[:current_spec_x],
        state[:current_spec_y],
        state[:current_spec_z],
        levels=g[:contourlevels],
        color=:grey50)
    g[:pltfit] = contour!(g[:axcontour],
        state[:current_fit_x],
        state[:current_fit_y],
        state[:current_fit_z],
        levels=g[:contourlevels],
        color=:orangered)

    # # create 3D plot
    # g[:ax3d] = Axis3(g[:panel3d][1,1], xlabel="δX / ppm", ylabel="δy / ppm", zlabel="Intensity",
    #     xreversed=true, yreversed=true)
    # g[:pltwireobs] = wireframe!(g[:ax3d],
    #     state[:current_spec_x],
    #     state[:current_spec_y],
    #     state[:current_spec_z],
    #     color=:grey80)
    # g[:pltwirefit] = wireframe!(g[:ax3d],
    #     state[:current_fit_x],
    #     state[:current_fit_y],
    #     state[:current_fit_z],
    #     color=:orangered)

    # add the peak positions
    g[:pltpeaks] = scatter!(g[:axcontour], state[:positions], markersize=10, marker=:x, color=state[:peakcolours])
    g[:pltlabels] = text!(g[:axcontour], state[:positions], text=state[:labels],
        fontsize=14,
        font=:bold,
        offset=(8,0),
        align=(:left,:center),
        color=:black)
    g[:pltinitialpeaks] = scatter!(g[:axcontour], state[:initialpositions], markersize=15, color=state[:peakcolours])

    # peak info
    g[:cmddelete] = Button(g[:panelinfo][1,1], label="Delete peak")
    g[:infotext] = Label(g[:panelinfo][2,1], state[:current_peak_info])

    # peak plot panel
    makepeakplot!(g, state, expt)

    addhanders!(g, state, expt)

    g[:fig]
end


function addhanders!(g, state, expt::FixedPeakExperiment)
    # link 2D and 3D contour plot limits
    if haskey(g, :ax3d)
        on(g[:axcontour].xaxis.attributes.limits) do xl
            xlims!(g[:ax3d], xl)
        end
        on(g[:axcontour].yaxis.attributes.limits) do yl
            ylims!(g[:ax3d], yl)
        end
    end

    # contour levels
    on(g[:cmdcontourup].clicks) do _
        g[:basecontour][] *= g[:contourscale][]
    end
    on(g[:cmdcontourdown].clicks) do _
        g[:basecontour][] /= g[:contourscale][]
    end

    # slices
    connect!(state[:current_slice], g[:sliderslice].value)
    on(g[:cmdsliceleft].clicks) do _
        i = state[:current_slice][]
        if i > 1
            i -= 1
            set_close_to!(g[:sliderslice], i)
        end
    end
    on(g[:cmdsliceright].clicks) do _
        i = state[:current_slice][]
        if i < nslices(expt)
            i += 1
            set_close_to!(g[:sliderslice], i)
        end
    end

    # peak hover
    onpick(g[:axcontour], g[:pltinitialpeaks]) do _, idx
        if state[:current_peak_idx][] != idx
            @debug "Setting current peak to $idx"
            state[:current_peak_idx][] = idx
            if haskey(g, :axpeakplot)
                autolimits!(g[:axpeakplot])
            end
        end
    end

    # delete peak
    on(g[:cmddelete].clicks) do _
        idx = state[:current_peak_idx][]
        if idx > 0
            state[:current_peak_idx][] = 0
            deletepeak!(expt, idx)
        end
    end

    # keyboard
    on(events(g[:axcontour]).keyboardbutton) do event
        if event.action == Keyboard.press && ispressed(g[:fig], Keyboard.a)
            @show pos = mouseposition(g[:axcontour])
            state[:total_peaks][] += 1
            id = "X$(state[:total_peaks][])"
            addpeak!(expt, Point2f(pos), id)
        end
    end
end


makepeakplot!(panel, state, expt::Experiment) = nothing