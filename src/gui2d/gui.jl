function gui!(expt::FixedPeakExperiment)
    GLMakie.activate!(; title="NMRAnalysis.jl (v$(string(pkgversion(GUI2D))))",
                      focus_on_show=true)

    state = expt.state[]

    g = Dict{Symbol,Any}() # GUI state
    state[:gui] = g

    g[:fig] = Figure(; size=(1200, 800))
    g[:paneltop] = g[:fig][1, 1:2] = GridLayout()
    g[:panelcontour] = g[:fig][2:3, 1] = GridLayout()
    g[:panelinfo] = g[:fig][2, 2] = GridLayout()
    g[:panelpeakplot] = g[:fig][3, 2] = GridLayout()
    rowsize!(g[:fig].layout, 1, Auto(true))
    rowsize!(g[:fig].layout, 2, Auto(false, 1))
    rowsize!(g[:fig].layout, 3, Auto(false, 1))
    colsize!(g[:fig].layout, 1, Auto(false, 2))
    colsize!(g[:fig].layout, 2, Auto(false, 1))

    # top panel
    g[:cmdcontourup] = Button(g[:paneltop][1, 1]; label="contour ↑")
    g[:cmdcontourdown] = Button(g[:paneltop][1, 2]; label="contour ↓")
    g[:cmdresetzoom] = Button(g[:paneltop][1, 3]; label="reset zoom")
    g[:sliderslice] = Slider(g[:paneltop][1, 4]; range=1:nslices(expt))
    g[:cmdsliceleft] = Button(g[:paneltop][1, 5]; label="←")
    g[:cmdsliceright] = Button(g[:paneltop][1, 6]; label="→")
    g[:slicelabel] = Label(g[:paneltop][1, 7], state[:current_slice_label])
    g[:togglefit] = Toggle(g[:paneltop][1, 8]; active=true)
    Label(g[:paneltop][1, 9], "Fitting")
    g[:cmdload] = Button(g[:paneltop][1, 10]; label="Load peak list")
    g[:cmdsave] = Button(g[:paneltop][1, 11]; label="Save to folder")

    # create contour plot
    g[:basecontour] = Observable(10.0)
    g[:contourscale] = Observable(1.7)
    g[:logscale] = lift(r -> [r .^ (0:10); -1 * (r .^ (0:10))], g[:contourscale])
    g[:contourlevels] = lift((c0, arr) -> c0 * arr, g[:basecontour], g[:logscale])

    g[:axcontour] = Axis(g[:panelcontour][1, 1];
                         xlabel="$(label(expt.specdata.nmrdata[1],F1Dim)) chemical shift (ppm)",
                         ylabel="$(label(expt.specdata.nmrdata[1],F2Dim)) chemical shift (ppm)",
                         xreversed=true, yreversed=true)
    g[:pltmask] = heatmap!(g[:axcontour],
                           state[:current_mask_x],
                           state[:current_mask_y],
                           state[:current_mask_z];
                           colormap=[:white, :lightgoldenrod1],
                           colorrange=(0, 1))
    g[:pltcontour] = contour!(g[:axcontour],
                              state[:current_spec_x],
                              state[:current_spec_y],
                              state[:current_spec_z];
                              levels=g[:contourlevels],
                              color=bicolours(:grey50, :lightblue))
    g[:pltfit] = contour!(g[:axcontour],
                          state[:current_fit_x],
                          state[:current_fit_y],
                          state[:current_fit_z];
                          levels=g[:contourlevels],
                          color=bicolours(:orangered, :dodgerblue))

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
    g[:pltpeaks] = scatter!(g[:axcontour], state[:positions]; markersize=10, marker=:x,
                            color=state[:peakcolours])
    g[:pltlabels] = text!(g[:axcontour], state[:positions]; text=state[:labels],
                          fontsize=14,
                          font=:bold,
                          offset=(8, 0),
                          align=(:left, :center),
                          color=:black)
    g[:pltinitialpeaks] = scatter!(g[:axcontour], state[:initialpositions]; markersize=15,
                                   color=state[:peakcolours])

    # peak info
    g[:cmdrename] = Button(g[:panelinfo][1, 1]; label="(R)ename peak")
    g[:cmddelete] = Button(g[:panelinfo][1, 2]; label="(D)elete peak")
    Label(g[:panelinfo][2, 1:2], "Press (A) to add new peak under mouse cursor";
          word_wrap=true)
    g[:sgradii] = SliderGrid(g[:panelinfo][3, 1:2],
                             (label="X radius", range=0.02:0.005:0.1, format="{:.3f} ppm",
                              startvalue=expt.xradius[]),
                             (label="Y radius", range=0.1:0.02:0.8, format="{:.2f} ppm",
                              startvalue=expt.yradius[])) # width = 350, tellheight = false)
    g[:sliderxradius] = g[:sgradii].sliders[1].value
    g[:slideryradius] = g[:sgradii].sliders[2].value

    g[:infotext] = Label(g[:panelinfo][4, 1:2], state[:current_peak_info])

    # peak plot panel
    makepeakplot!(g, state, expt)

    @debug "Adding handlers"
    addhanders!(g, state, expt)

    @debug "Showing figure"
    return g[:fig]
end

function addhanders!(g, state, expt::FixedPeakExperiment)
    g[:fig].scene.backgroundcolor = lift(state[:mode]) do mode
        if mode == :fitting
            RGBAf(1.0, 0.63, 0.48, 1.0)      # :salmon
        elseif mode == :renaming || mode == :renamingstart
            RGBAf(0.75, 0.94, 1.0, 1.0)     # :lightblue
        elseif mode == :moving
            RGBAf(0.6, 0.98, 0.6, 1.0)      # :palegreen
        else
            RGBAf(1.0, 1.0, 1.0, 1.0)       # :white
        end
    end

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

    # zoom
    on(g[:cmdresetzoom].clicks) do _
        reset_limits!(g[:axcontour])
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

    # fitting active
    connect!(expt.isfitting, g[:togglefit].active)
    connect!(g[:pltfit].visible, g[:togglefit].active)

    # load peak list
    on(g[:cmdload].clicks) do _
        loadpeaks!(expt)
    end

    # save peak list
    on(g[:cmdsave].clicks) do _
        saveresults!(expt)
    end

    # delete peak
    on(g[:cmddelete].clicks) do _
        idx = state[:current_peak_idx][]
        if idx > 0
            state[:current_peak_idx][] = 0
            deletepeak!(expt, idx)
        end
    end

    # rename peak
    on(g[:cmdrename].clicks) do _
        if state[:current_peak_idx][] > 0
            renamepeak!(expt, state, :mouse)
        end
    end

    # radius sliders
    connect!(expt.xradius, g[:sliderxradius])
    connect!(expt.yradius, g[:slideryradius])

    # peak hover
    onpick(g[:axcontour], g[:pltinitialpeaks]) do _, idx
        # if state[:current_peak_idx][] != idx && state[:mode][] == :normal
        if state[:mode][] == :normal
            @debug "Setting current peak to $idx"
            state[:current_peak_idx][] = idx
            if haskey(g, :axpeakplot)
                autolimits!(g[:axpeakplot])
            end
            if haskey(g, :axpeakplotX)
                autolimits!(g[:axpeakplotX])
            end
            if haskey(g, :axpeakplotY)
                autolimits!(g[:axpeakplotY])
            end
        end
    end
    # mouse handling
    on(events(g[:axcontour]).mousebutton; priority=2) do event
        process_mousebutton(expt, state, event)
    end
    on(events(g[:axcontour]).mouseposition; priority=2) do mousepos
        process_mouseposition(expt, state, mousepos)
    end

    # keyboard
    on(events(g[:axcontour]).keyboardbutton; priority=2) do event
        process_keyboardbutton(expt, state, event)
    end
    on(events(g[:fig]).unicode_input) do character
        process_unicode_input(expt, state, character)
    end
end

function renamepeak!(expt, state, initiator)
    @debug "Renaming peak"
    if initiator == :keyboard
        state[:mode][] = :renamingstart
    else
        state[:mode][] = :renaming
    end
    # state[:current_peak_info][] = "Renaming peak $(peak.label[]), press Enter to finish, Backspace to delete, Esc to cancel"
    state[:oldlabel][] = state[:current_peak][].label[]
    state[:current_peak][].label[] = "‸"
    notify(expt.peaks)
end

bicolours(c1, c2) = [fill(c1, 11); fill(c2, 11)]