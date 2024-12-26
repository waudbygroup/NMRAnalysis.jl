function loadpeaklist(peaklistfilename, spectra)
    return Observable(Vector{Peak2D}())
end

struct Setup2D
    peakradiusX
    peakradiusY

    maxpeakmovementX
    maxpeakmovementY

    minR2
    maxR2
end

Setup2D() = Setup2D(Observable(0.03), Observable(0.15), 0.03, 0.15, 2., 200.)
struct ControlPanel
    □contourup
    □contourdown
    □fitting
    □expfitting
    □outputfilename
    □savelist
    □info
    □plane
    □Xradius
    □Yradius
end
function ControlPanel(fig)
    □planeslider = labelslider!(fig, "plane: ", 1:1, startvalue=1)

    □contourup = Button(fig, label="contours ↑")
    □contourdown = Button(fig, label="contours ↓")

    □Xradiusslider = labelslider!(fig, "X fitting radius: ", 0.005:0.005:0.1, startvalue=0.05)
    □Yradiusslider = labelslider!(fig, "Y fitting radius: ", 0.05:0.05:1, startvalue=0.5)

    □fitting = Toggle(fig, active=false)
    □expfitting = Toggle(fig, active=false)
    □outputfilename = Textbox(fig, placeholder="Enter output filename", stored_string="peak-fits.txt")#, validator=String)
    □savelist = Button(fig, label="save peak list")
    □info = Label(fig, "", justification=:center, lineheight=1.2)
    
    return ControlPanel(□contourup, □contourdown, □fitting, □expfitting, □outputfilename, □savelist, □info, □planeslider, □Xradiusslider, □Yradiusslider)
end

struct PlotPanel
    axis
    contourplot
    scatterplot
    labelsplot
    maskplot
    fitplot
end
PlotPanel() = PlotPanel(nothing, nothing, nothing, nothing, nothing, nothing)

struct FitPanel
    axis
    seriesplot
    fitplot
end
FitPanel() = FitPanel(nothing, nothing, nothing)


function fit2d(spectra, fit_function=nothing, peaklistfilename=nothing)
    # parse input peak list
    peaks = loadpeaklist(peaklistfilename, spectra)

    specdata = preparespecdata(spectra)

    # setup parameters
    setup = Setup2D()

    # create state
    fig = Figure()
    state = State(fig, setup, spectra, peaks)

    # create figure
    preparefigure!(state)

    # set up callbacks
    # - contour levels
    on(state.controls.□contourup.clicks) do clicks
        state.baselev[] *= 1.3
    end
    on(state.controls.□contourdown.clicks) do clicks
        state.baselev[] /= 1.3
    end
    # - mouse dragging
    on(events(state.fig).mousebutton, priority = 2) do event
        process_mousebutton(state, event)
    end
    on(events(state.fig).mouseposition, priority = 2) do mousepos
        process_mouseposition(state, mousepos)
    end
    
    # - key commands
    on(events(state.fig).keyboardbutton, priority = 2) do event
        process_keyboardbutton(state, event)
    end
    on(events(state.fig).unicode_input) do character
        process_unicode_input(state, character)
    end
    # # - fitting
    onany(state.peaks, state.controls.□fitting.active) do P, fitting
        if !fitting
            state.controls.□expfitting.active[] = false
        end
        fitting || return Consume(false)
        state.dragging[] && return Consume(false)

        fit!(state)
        return Consume()
    end
    # # - exp fitting
    onany(state.peaks, state.controls.□expfitting.active) do P, expfitting
        expfitting || return Consume(false)
        state.dragging[] && return Consume(false)

        state.controls.□expfitting.active[] || return Consume(false)
        fitexp!(state)

        return Consume()
    end
    # - save peak list
    on(state.controls.□savelist.clicks) do clicks
        savepeaklist(state)
    end

    # display app
    return state.fig
end

function preparefigure!(state)
    mainpanel = state.fig[1,1]
    state.fig[1, 1] = GridLayout()
    controllayout = state.fig[1, 2] = GridLayout()
    
    preparecontrolpanel!(state, controllayout)
    preparecontourpanel!(state, mainpanel)
    preparefitplot!(state, fitplotpanel)

    # only display fitted spectrum if fitting is enabled
    connect!(state.plots.fitplot.alpha, state.controls.□fitting.active)

    # only display exponential fits if exponential fitting is enabled
    connect!(state.fitplots.fitplot.visible, state.controls.□expfitting.active)

    bgcolor = lift(x -> x ? colorant"pink" : colorant"white", state.renaming)
    connect!(state.fig.scene.backgroundcolor, bgcolor)
    
    infolabeltext = getpeakinfo(state)
    connect!(state.controls.□info.text, infolabeltext)

    connect!(state.currentslice, state.controls.□plane.slider.value)
    connect!(state.setup.peakradiusX, state.controls.□Xradius.slider.value)
    connect!(state.setup.peakradiusY, state.controls.□Yradius.slider.value)
end


function preparecontourpanel!(state, panel)
    clev = state.baselev[] * exp.(range(0,5,11))
    contourax = Axis(panel,
        xreversed=true, yreversed=true,
        xlabel=lift(i->label(state.spectra[i], F1Dim), state.currentslice),
        ylabel=lift(i->label(state.spectra[i], F2Dim), state.currentslice),
        title=lift(i->choptitle(label(state.spectra[i])), state.currentslice),
        xgridvisible=false, ygridvisible=false)
    contourplot = lift(state.currentslice) do i
        x = state.X[i]
        y = state.Y[i]
        z = state.Z[i]
        
        contour!(contourax, x, y, z, levels=clev,
            colormap=[colorant"red", colorant"black"],
            colorrange=(-0.001*state.specdata.σ, 0.001*state.specdata.σ),
            lowclip=:red,
            highclip=:black,
            inspectable=false)
    end
    
    # currentZsim = lift(state.specdata.Zfit, state.controls.□plane.slider.value) do Zdata, slice
    #     Zdata[slice]
    # end

    

    maskplot = heatmap!(contourax, state.specdata.dX, state.specdata.dY, getmask(state),
        colormap=[colorant"white", colorant"navajowhite"], inspectable=false)
    fitplot = contour!(contourax,
        state.specdata.dX,
        state.specdata.dY,
        currentZsim,
        levels=state.specdata.clev,
        colormap=[colorant"cyan", colorant"magenta"],
        colorrange=(-0.001*state.specdata.σ, 0.001*state.specdata.σ),
        lowclip=:cyan,
        highclip=:magenta,
        inspectable=false)
    scatterplot = scatter!(contourax, getpeakpositions(state),
        color=getpeakcolors(state),
        markersize=15,
        marker=:circle,
        inspectable=true)
    labelsplot = text!(contourax, getpeaktext(state), offset=(10,10))

    inspector = DataInspector(contourax)

    state.plots = PlotPanel(contourax, contourplot)#, scatterplot, labelsplot, maskplot, fitplot)
end


function preparecontrolpanel!(state, layout)
    layout[1,1] = state.controls.□plane.layout
    state.controls.□plane.slider.range[] = 1:state.specdata.N
    layout[2,1] = grid!([state.controls.□contourup ;; state.controls.□contourdown])
    layout[3,1] = state.controls.□Xradius.layout
    layout[4,1] = state.controls.□Yradius.layout
    layout[5,1] = grid!([Label(state.fig, "peak fitting") ;; state.controls.□fitting])
    layout[6,1] = grid!([Label(state.fig, "exponential fitting") ;; state.controls.□expfitting])
    layout[7,1] = grid!([state.controls.□savelist ;; state.controls.□outputfilename])
    layout[8,1] = state.controls.□info
    
    layout.tellheight[] = true
end