function integrategui!(state)
    gui = state["gui"]
    gui["shouldclose"] = false

    fig = Figure(size=(1200, 800))

    # spectrum display
    spectraax = Axis(fig[1, 1], xlabel="Chemical shift (ppm)", ylabel="Intensity",
        xreversed=true,
        yzoomlock=true,
        ypanlock=true,
        xrectzoom=true,
        yrectzoom=false,
        xgridvisible=false,
        ygridvisible=false,)
    hideydecorations!(spectraax)
    vspan!(spectraax, state["standardx1"], state["standardx2"], color=Makie.wong_colors()[1], alpha=0.3)
    vspan!(spectraax, state["unknownx1"], state["unknownx2"], color=Makie.wong_colors()[2], alpha=0.3)
    series!(spectraax, gui["linedata"], color=gui["linecolors"])

    # side panel
    side_panel = fig[1, 2] = GridLayout()

    top_row = side_panel[1, 1] = GridLayout()
    button_up = top_row[1, 1] = Button(fig, label = "↑")
    button_down = top_row[1, 2] = Button(fig, label = "↓")
    button_zoomout = top_row[1, 3] = Button(fig, label = "Reset zoom")
    button_quit = top_row[1, 5] = Button(fig, label = "Quit")

    second_row = side_panel[2, 1] = GridLayout()
    text_spectra_filename = second_row[1, 1] = Textbox(fig, stored_string = gui["spectra_filename"], width=150, tellwidth=true)
    button_savespectra = second_row[1, 2] = Button(fig, label = "Save spectra plot")

    standard_row = side_panel[3, 1] = GridLayout()
    standard_row[1, 1] = Label(fig, gui["integration_label_1"])
    button_setstandard = standard_row[1, 2] = Button(fig, label = "Set to display region")

    unknown_row = side_panel[4, 1] = GridLayout()
    unknown_row[1, 1] = Label(fig, gui["integration_label_2"])
    button_setunknown = unknown_row[1, 2] = Button(fig, label = "Set to display region")
    
    save_row = side_panel[5, 1] = GridLayout()
    text_results_filename = save_row[1, 1] = Textbox(fig, stored_string = gui["results_filename"], width=150, tellwidth=true)
    button_saveresults = save_row[1, 2] = Button(fig, label = "Save results")

    blank_row = side_panel[6, 1] = GridLayout()

    heatmap_row = side_panel[7, 1] = GridLayout()
    heatmap_ax = Axis(heatmap_row[1,1], xlabel="", ylabel="", yreversed=true, xticks=1:12, yticks=(1:8, ["A", "B", "C", "D", "E", "F", "G", "H"]))
    hm = heatmap!(heatmap_ax, gui["integralplate"], colormap=:viridis)
    scatter!(heatmap_ax, gui["selected_points"], markersize=10, color=gui["linecolors"])
    Colorbar(heatmap_row[1,2], hm, label="Integral ratio")

    side_panel[8, 1] = Label(fig, """
SPACE - toggle selection; 'a' - select all
'r' - select row; 'c' - select column
SHIFT + 'r' - add row; SHIFT + 'c' - add column
    """)

    # event handling
    on(button_up.clicks) do _
        gui["plotscale"][] *= 1.3
        @show gui["plotscale"][]
    end
    on(button_down.clicks) do _
        gui["plotscale"][] /= 1.3
    end
    on(button_zoomout.clicks) do _
        gui["plotscale"][] = 1
        autolimits!(spectraax)
    end
    on(text_spectra_filename.stored_string) do s
        gui["spectra_filename"][] = s
    end
    on(text_results_filename.stored_string) do s
        gui["results_filename"][] = s
    end
    on(button_quit.clicks) do _
        gui["shouldclose"] = true
    end
    on(button_setstandard.clicks) do _
        xl1, xl2 = spectraax.xaxis.attributes.limits[]
        x = (xl1 + xl2) / 2
        dx = (xl2 - xl1)
        state["standardx"][] = x
        state["standarddx"][] = dx
    end
    on(button_setunknown.clicks) do _
        xl1, xl2 = spectraax.xaxis.attributes.limits[]
        x = (xl1 + xl2) / 2
        dx = (xl2 - xl1)
        state["unknownx"][] = x
        state["unknowndx"][] = dx
    end
    on(button_savespectra.clicks) do _
        filename = joinpath(state["foldername"], gui["spectra_filename"][])
        # save(filename, spectraax.scene; backend=CairoMakie)
        save(filename, fig; backend=CairoMakie)
    end
    on(button_saveresults.clicks) do _
        saveresults(state)
    end

    # keyboard
    on(events(fig.scene).keyboardbutton) do event
        # check if mouse is within heatmap (for later)
        x, y = mouseposition(heatmap_ax)
        i, j = round(Int, x), round(Int, y)
        inheatmap = (1 <= i <= 12 && 1 <= j <= 8)

        if event.action == Keyboard.press || event.action == Keyboard.repeat
            if event.key == Keyboard.up
                gui["plotscale"][] *= 1.3
            elseif event.key == Keyboard.down
                gui["plotscale"][] /= 1.3
            elseif event.key == Keyboard.left
                xl1, xl2 = spectraax.xaxis.attributes.limits[]
                dx = (xl2 - xl1) * 0.1
                xlims!(spectraax, xl2 - dx, xl1 - dx)
            elseif event.key == Keyboard.right
                xl1, xl2 = spectraax.xaxis.attributes.limits[]
                dx = (xl2 - xl1) * 0.1
                xlims!(spectraax, xl2 + dx, xl1 + dx)
            end
        end
        if inheatmap && event.action == Keyboard.press
            if ispressed(fig, Keyboard.space)
                # invert selection
                sel = state["selected"][][i + 12*(j-1)]
                if sel && state["nselected"][] > 1 # check - don't switch off all spectra
                    state["selected"][][i + 12*(j-1)] = false
                    notify(state["selected"])
                    autolimitsy!(spectraax)
                end
                if !sel
                    state["selected"][][i + 12*(j-1)] = true
                    notify(state["selected"])
                    autolimitsy!(spectraax)
                end
            elseif ispressed(fig, Keyboard.a)
                # select all
                for col=1:12, row=1:8
                    state["selected"][][col + 12*(row-1)] = true
                end
                notify(state["selected"])
                autolimitsy!(spectraax)
            elseif ispressed(fig, Exclusively(Keyboard.r))
                # select row
                for col=1:12, row=1:8
                    state["selected"][][col + 12*(row-1)] = row == j
                end
                notify(state["selected"])
                autolimitsy!(spectraax)
            elseif ispressed(fig, Exclusively(Keyboard.c))
                # select column
                for col=1:12, row=1:8
                    state["selected"][][col + 12*(row-1)] = col == i
                end
                notify(state["selected"])
                autolimitsy!(spectraax)
            elseif ispressed(fig, Keyboard.r & (Keyboard.left_shift | Keyboard.right_shift))
                # add row
                for col=1:12
                    state["selected"][][col + 12*(j-1)] = true
                end
                notify(state["selected"])
                autolimitsy!(spectraax)
            elseif ispressed(fig, Keyboard.c & (Keyboard.left_shift | Keyboard.right_shift))
                # add column
                for row=1:8
                    state["selected"][][i + 12*(row-1)] = true
                end
                notify(state["selected"])
                autolimitsy!(spectraax)
            end
        end
    end

    # mouse
    on(events(spectraax).scroll) do (_, dy)
        if dy < 0
            gui["plotscale"][] *= 1.3
        else
            gui["plotscale"][] /= 1.3
        end
    end
    

    display(fig)
    while !gui["shouldclose"]
        sleep(0.1)
        if !isopen(fig.scene)
            break
        end
    end
    
    GLMakie.closeall()
end


function autolimitsy!(ax)
    xl1, xl2 = ax.xaxis.attributes.limits[]
    autolimits!(ax)
    xlims!(ax, xl2, xl1)
end