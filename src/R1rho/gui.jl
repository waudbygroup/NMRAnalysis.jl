function gui!(state)
    state[:gui] = Dict{Symbol, Any}()
    gui = state[:gui]
    fig = Figure(size=(1200,800))

    c1 = Makie.wong_colors()[1]
    c2 = Makie.wong_colors()[2]
    c3 = Makie.wong_colors()[3]
    c4 = Makie.wong_colors()[4]
    c5 = Makie.wong_colors()[5]
    c6 = Makie.wong_colors()[6]
    c7 = Makie.wong_colors()[7]

    top_panel = fig[1, 1] = GridLayout()
    bottom_panel = fig[2, 1] = GridLayout()

    gui[:specplottitle] = lift(state[:currentseries]) do i
        "Observed spectrum (νSL = $(round(0.001*νSL(state[:dataset])[i],digits=2)) kHz)"
    end
    ax_spectra = Axis(top_panel[1,1],
        xreversed=true,
        yzoomlock=true,
        ypanlock=true,
        xrectzoom=true,
        yrectzoom=false,
        xlabel="Chemical shift (ppm)",
        ylabel="Intensity",
        title=gui[:specplottitle]
        )
    hideydecorations!(ax_spectra)

    gui[:noisespan] = lift(state[:noiseppm], state[:dx]) do noise, dx
        noise-dx/2, noise+dx/2
    end
    gui[:peakspan] = lift(state[:peakppm], state[:dx]) do peak, dx
        peak-dx/2, peak+dx/2
    end
    hlines!(ax_spectra, [0], color=:grey)
    gui[:peakspan] = vspan!(ax_spectra, @lift($(gui[:peakspan])[1]), @lift($(gui[:peakspan])[2]), alpha=0.7, label="Peak")
    gui[:noisespan] = vspan!(ax_spectra, @lift($(gui[:noisespan])[1]), @lift($(gui[:noisespan])[2]), alpha=0.7, label="Noise", color=c4)
    lines!(ax_spectra, state[:currentspectrum], label = "Observed")
    axislegend(ax_spectra, position=:lt)

    input_panel = top_panel[1,2] = GridLayout()

    input_panel[1,1] = Label(fig, "Series:")
    slider_current = input_panel[1,2] = Slider(fig, range = 1:state[:nseries], width=150)
    gui[:slider_current] = slider_current
    connect!(state[:currentseries], slider_current.value)

    input_panel[2,1] = Label(fig, "Peak position (ppm):")
    text_peakppm = input_panel[2,2] = Textbox(fig, stored_string=string(round(state[:peakppm][],digits=2)), validator=Float64, width=150)
    gui[:text_peakppm] = text_peakppm
    on(text_peakppm.stored_string) do s
        state[:peakppm][] = parse(Float64, s)
    end

    input_panel[3,1] = Label(fig, "Noise position (ppm):")
    text_noiseppm = input_panel[3,2] = Textbox(fig, stored_string=string(round(state[:noiseppm][],digits=2)), validator=Float64, width=150)
    gui[:text_noiseppm] = text_noiseppm
    on(text_noiseppm.stored_string) do s
        state[:noiseppm][] = parse(Float64, s)
    end

    input_panel[4,1] = Label(fig, "Integration width (ppm):")
    text_dx = input_panel[4,2] = Textbox(fig, stored_string=string(round(state[:dx][], digits=2)), validator=Float64, width=150)
    gui[:text_dx] = text_dx
    on(text_dx.stored_string) do s
        state[:dx][] = parse(Float64, s)
    end

    # input_panel[5,1] = Label(fig, "Fit")
    # cb_fit = input_panel[5, 2] = Checkbox(fig, checked = state[:isfitting][])
    # connect!(state[:isfitting], cb_fit.checked)

    input_panel[5,1] = Label(fig, "Initial I0:")
    text_I0 = input_panel[5,2] = Textbox(fig, stored_string=string(round(state[:initialI0][], digits=1)), validator=Float64, width=150)
    gui[:text_I0] = text_I0
    on(text_I0.stored_string) do s
        state[:initialI0][] = parse(Float64, s)
    end
    on(state[:intensities]) do I
        state[:initialI0][] = I[1]
        text_I0.displayed_string[] = string(round(I[1], digits=1))
    end

    input_panel[6,1] = Label(fig, "Initial R2,0 (s⁻¹):")
    text_R20 = input_panel[6,2] = Textbox(fig, stored_string=string(round(state[:initialR20][], digits=1)), validator=Float64, width=150)
    gui[:text_R20] = text_R20
    on(text_R20.stored_string) do s
        state[:initialR20][] = parse(Float64, s)
    end

    input_panel[7,1] = Label(fig, "Initial Rex (s⁻¹):")
    text_Rex = input_panel[7,2] = Textbox(fig, stored_string=string(round(state[:initialRex][], digits=1)), validator=Float64, width=150)
    gui[:text_Rex] = text_Rex
    on(text_Rex.stored_string) do s
        state[:initialRex][] = parse(Float64, s)
    end

    input_panel[8,1] = Label(fig, "Initial kex (s⁻¹):")
    text_kex = input_panel[8,2] = Textbox(fig, stored_string=string(round(exp(state[:initiallnk][]), digits=1)), validator=Float64, width=150)
    gui[:text_kex] = text_kex
    on(text_kex.stored_string) do s
        state[:initiallnk][] = log(parse(Float64, s))
    end

    gui[:fitplottitle] = lift(state[:currentseries]) do i
        "Peak integrals (νSL = $(round(0.001*νSL(state[:dataset])[i],digits=2)) kHz)"
    end
    ax_fit = Axis(bottom_panel[1,1],
        xlabel = "TSL (ms)",
        ylabel = "Peak integral",
        title = gui[:fitplottitle]
        )
    errorbars!(ax_fit, state[:currenterror])
    scatter!(ax_fit, state[:currentscatter], label = "Observed")
    lines!(ax_fit, state[:currentfit], label = "Global fit")
    lines!(ax_fit, state[:currentexpfit], label = "Exponential fit")
    axislegend(ax_fit, position=:rt)

    ax_fit_R1rho = Axis(bottom_panel[1,2],
        xlabel = "νSL (kHz)",
        ylabel = "R1rho (s⁻¹)",
        title = "Dispersion curve"
        )
    errorbars!(ax_fit_R1rho, state[:expfiterror], color=c2)
    scatter!(ax_fit_R1rho, state[:expfitpoints], label = "Exponential fits", color=c2)
    lines!(ax_fit_R1rho, state[:fitR1rho], label = "Global fit")
    axislegend(ax_fit_R1rho, position=:rt)

    results_panel = bottom_panel[1,3] = GridLayout(tellheight=false)
    Label(results_panel[1,1:2], lift(x->"I0: $x", state[:fitI0]))
    Label(results_panel[2,1:2], lift(x->"R2,0 (s⁻¹): $x", state[:fitR20]))
    Label(results_panel[3,1:2], lift(x->"Rex (s⁻¹): $x", state[:fitRex]))
    Label(results_panel[4,1:2], lift(x->"kex (s⁻¹): $(exp(x))", state[:fitlnk]))
    results_panel[5,1:2] = Label(fig, "Working directory:\n$(pwd())")
    results_panel[6,1] = Label(fig, "Output folder:")
    text_out = results_panel[6,2] = Textbox(fig, stored_string="out", width=150)
    gui[:text_out] = text_out
    on(text_out.stored_string) do s
        state[:outputdir][] = s
    end
    button_save = results_panel[7,1:2] = Button(fig, label="Save results")
    on(button_save.clicks) do _
        savefig!(state)
    end

    dragging = :nothing
    on(events(ax_spectra).mousebutton) do event
        global dragging, dragidx
        if event.button == Mouse.left
            if event.action == Mouse.press
                if mouseover(fig, gui[:peakspan])
                    dragging = :peak
                elseif mouseover(fig, gui[:noisespan])
                    dragging = :noise
                else
                    dragging = :nothing
                end
                return Consume(dragging != :nothing)
            elseif event.action == Mouse.release
                # Exit dragging
                dragging = :nothing
                return Consume(false)
            end    
        end
    end
    on(events(fig).mouseposition, priority = 2) do mp
        global dragging
        if dragging != :nothing
            p = mouseposition(ax_spectra)
            if dragging == :peak
                state[:peakppm][] = p[1]
            elseif dragging == :noise
                state[:noiseppm][] = p[1]
            end
            return Consume(true)
        end
        return Consume(false)
    end
    

    display(fig)
    # while !state["should_close"][]
    #     sleep(0.1)
    #     if !isopen(fig.scene)
    #         break
    #     end
    # end
    
    # GLMakie.closeall()
end



function savefig!(state)
    outputdir = state[:outputdir][]
    @info "Saving results to $(joinpath(pwd(),outputdir))"
    if !isdir(outputdir)
        mkdir(outputdir)
    end

    c1 = Makie.wong_colors()[1]
    c2 = Makie.wong_colors()[2]
    
    # dispersion fit
    fig = Figure()
    ax_fit_R1rho = Axis(fig[1,1],
        xlabel = "νSL (kHz)",
        ylabel = "R1rho (s⁻¹)"
        )
    errorbars!(ax_fit_R1rho, state[:expfiterror], color=c2)
    scatter!(ax_fit_R1rho, state[:expfitpoints], label = "Exponential fits", color=c2)
    lines!(ax_fit_R1rho, state[:fitR1rho], label = "Global fit")
    axislegend(ax_fit_R1rho, position=:rt)
    save(joinpath(outputdir, "dispersion.pdf"), fig; backend=CairoMakie)

    # intensities
    for i=1:state[:nseries]
        title = "Peak integrals (νSL = $(round(0.001*νSL(state[:dataset])[i],digits=2)) kHz)"
        filename = "intensities_$(round(0.001*νSL(state[:dataset])[i],digits=2))_kHz.pdf"
        fig = Figure()
        ax_fit = Axis(fig[1,1],
            xlabel = "TSL (ms)",
            ylabel = "Peak integral",
            title = title
            )
        errorbars!(ax_fit, state[:errorpoints][][i])
        scatter!(ax_fit, state[:scatterpoints][][i], label = "Observed")
        lines!(ax_fit, state[:fitseries][][i], label = "Global fit")
        axislegend(ax_fit, position=:rt)
        save(joinpath(outputdir, filename), fig; backend=CairoMakie)
    end

    # save fit results to CSV: peak/noise positions, and fit parameters. don't use additional libraries
    filename = joinpath(outputdir, "results.txt")
    open(filename, "w") do f
        for filename in state[:filenames]
            println(f, "Input file: $filename")
        end
        println(f, "Peak position (ppm):, $(state[:peakppm][])")
        println(f, "Noise position (ppm):, $(state[:noiseppm][])")
        println(f, "")
        println(f, "Integration width (ppm):, $(state[:dx][])")
        println(f, "Initial I0:, $(state[:initialI0][])")
        println(f, "Initial R2,0 (s⁻¹):, $(state[:initialR20][])")
        println(f, "Initial Rex (s⁻¹):, $(state[:initialRex][])")
        println(f, "Initial kex (s⁻¹):, $(exp(state[:initiallnk][]))")
        println(f, "")
        println(f, "Fitted I0:, $(state[:fitI0][])")
        println(f, "Fitted R2,0 (s⁻¹):, $(state[:fitR20][])")
        println(f, "Fitted Rex (s⁻¹):, $(state[:fitRex][])")
        println(f, "Fitted kex (s⁻¹):, $(exp(state[:fitlnk][]))")
    end
end

