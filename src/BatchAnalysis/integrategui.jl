function integrategui!(state)
    gui = state["gui"]

    fig = Figure(size=(800, 600))

    top_panel = fig[1, 1] = GridLayout()
    button_up = top_panel[1, 1] = Button(fig, label = "↑")
    button_down = top_panel[1, 2] = Button(fig, label = "↓")
    button_zoomout = top_panel[1, 3] = Button(fig, label = "Reset zoom")
    top_panel[1, 4] = Label(fig, "Ref. chemical shift (ppm):")
    text_refshift = top_panel[1, 5] = Textbox(fig, stored_string = string(state["referenceshift"][]),
        validator = Float64, tellwidth = false, width=50)
    button_setref = top_panel[1, 6] = Button(fig, label = "Set reference to visible range")
    button_cancel = top_panel[1, 7] = Button(fig, label = "Cancel")
    button_continue = top_panel[1, 8] = Button(fig, label = "Continue")

    spectraax = Axis(fig[2, 1], xlabel="Chemical shift (ppm)", ylabel="Intensity",
        xreversed=true,
        yzoomlock=true,
        ypanlock=true,
        xrectzoom=true,
        yrectzoom=false)
    hideydecorations!(spectraax)
    spectraplot = series!(spectraax, gui["linedata"], color=gui["linecolors"])

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
        reset_limits!(spectraax)
    end
    
    on(text_refshift.stored_string) do s
        state["referenceshift"][] = parse(Float64, s)
    end

    # keyboard
    on(events(fig.scene).keyboardbutton) do event
        if event.action == Keyboard.press || event.action == Keyboard.repeat
            if event.key == Keyboard.up
                gui["plotscale"][] *= 1.3
            elseif event.key == Keyboard.down
                gui["plotscale"][] /= 1.3
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
    

    wait(display(fig))
end