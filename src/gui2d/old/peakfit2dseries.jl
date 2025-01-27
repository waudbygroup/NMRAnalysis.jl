struct Setup
    peakradiusX
    peakradiusY

    maxpeakmovementX
    maxpeakmovementY

    minR2
    maxR2
end
Setup() = Setup(Observable(0.03), Observable(0.15), 0.03, 0.15, 2., 200.)

struct SpecData
    spectrum
    dX
    dY
    T
    Z
    σ

    Zfit

    clev0
    clev
end

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
    # □planeslider = labelslider!(fig, "plane: ", 1:1, startvalue=1)
    □planeslider = Slider(fig, range = 1:1, startvalue=1)

    □contourup = Button(fig, label="contours ↑")
    □contourdown = Button(fig, label="contours ↓")

    # □Xradiusslider = labelslider!(fig, "X fitting radius: ", 0.005:0.005:0.1, startvalue=0.05)
    # □Yradiusslider = labelslider!(fig, "Y fitting radius: ", 0.05:0.05:1, startvalue=0.5)
    □Xradiusslider = Slider(fig, range=0.005:0.005:0.1, startvalue=0.05)
    □Yradiusslider = Slider(fig, range=0.05:0.05:1, startvalue=0.5)

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

mutable struct State
    fig

    setup::Setup
    specdata::SpecData
    peaks::Observable{Vector{PeakPseudo2D}}

    clusters::Vector{Vector{Int64}}

    dragging::Observable{Bool}
    renaming::Observable{Bool}
    dragpeakindex::Int64
    activepeakindex::Int64
    highlightedpeakindex::Observable{Int64}
    oldlabel::String
    startinglabelinput::Bool

    plots::PlotPanel
    controls::ControlPanel
    fitplots::FitPanel
end
State(fig, setup, specdata, peaks) = State(fig, setup, specdata, peaks,
        makeclusters(peaks),
        Observable(false), Observable(false), 0, 0, Observable(0), "", false,
        PlotPanel(),
        ControlPanel(fig),
        FitPanel())





function loadspectrum(inputfilename, taufilename)
    spec = loadnmr(inputfilename)
    @show T = vec(readdlm(taufilename, comments=true))
    
    dX = data(spec, F1Dim)
    dY = data(spec, F2Dim)
    
    Z = data(spec) / scale(spec)
    σ = spec[:noise] / scale(spec)

    Zfit = Observable(0Z)

    # prepare contour levels
    baselev = exp.(LinRange(0,5,11))
    baselev = sort([baselev; -baselev])
    clev0 = Observable(5σ)
    clev = @lift $clev0 * baselev
    @debug "Initial contour levels" clev

    return SpecData(spec, dX, dY, T, Z, σ, Zfit, clev0, clev)
end



function loadpeaklist(peaklistfilename, specdata)
    isnothing(peaklistfilename) && return Observable(Vector{PeakPseudo2D}())

    nt = length(specdata.T)
    @show amp0 = zeros(nt)

    peaklist = readdlm(peaklistfilename, comments=true)

    npeaks, ncol = size(peaklist)
    ncol < 3 && return Observable(Vector{PeakPseudo2D}())

    peaks = Vector{PeakPseudo2D}()
    for i=1:npeaks
        label = peaklist[i,1]
        dX = peaklist[i,2]
        dY = peaklist[i,3]
        # if ncol > 4
        #     lwX = (ismissing(peaklist[i,4]) || peaklist[i,4]=="missing") ? 30. : peaklist[i,4]
        #     lwY = (ismissing(peaklist[i,5]) || peaklist[i,5]=="missing") ? 30. : peaklist[i,5]
        # end
        # if ncol > 5
        #     amp = (ismissing(peaklist[i,6]) || peaklist[i,6]=="missing") ? [] : peaklist[i,6:end]
        # end
        # push!(peaks, PeakPseudo2D(Point2f(dX, dY), amp, lwX, lwY, label))
        push!(peaks, PeakPseudo2D(Point2f(dX, dY), amp0, 30., 30., label))
        # PeakPseudo2D(newpos, amp, 30, 30, "X$(length(state.peaks[])+1)")
    end

    return Observable(peaks)
end



getpeakpositions(state) = lift(P->[p.position for p ∈ P], state.peaks)
getpeakcolors(state) = lift(P->[p.touched ? "orangered" : "dodgerblue" for p ∈ P], state.peaks)
function getpeaktext(state)
    lift(state.peaks) do P
        if isempty(P)
            [("",
             Point2f((state.specdata.dX[1]+state.specdata.dX[end])/2,
                     (state.specdata.dY[1]+state.specdata.dY[end])/2))]
        else
            [(p.label, p.position) for p ∈ P]
        end
    end
end
function getmask(state)
    peakpositions = getpeakpositions(state)
    lift(peakpositions, state.setup.peakradiusX, state.setup.peakradiusY) do peakpos, xw, yw
        m = zeros(Bool, size(state.specdata.Z, 1), size(state.specdata.Z, 2))
        for pos ∈ peakpos
            maskellipse!(m, state.specdata.dX, state.specdata.dY, pos[1], pos[2], xw, yw)
        end
        m
    end
end

function getpeakinfo(state)
    lift(state.highlightedpeakindex, state.peaks) do i, P
        if i==0 || i > length(P)
            ""
        else
            p = P[i]
            if p.touched
                """
                $(p.label)
                δ(X): $(round(p.position[1], sigdigits=4)) ppm
                δ(Y): $(round(p.position[2], sigdigits=4)) ppm
                """
            else
                if state.controls.□expfitting.active[]
                    """
                    $(p.label)
                    δ(X): $(p.position[1] ± p.position_err[1]) ppm
                    δ(Y): $(p.position[2] ± p.position_err[2]) ppm
                    lw(X): $(p.lwX ± p.lwX_err) Hz
                    lw(Y): $(p.lwY ± p.lwY_err) Hz
                    relaxation rate: $(p.k ± p.k_err) s⁻¹
                    """
                else
                    """
                    $(p.label)
                    δ(X): $(p.position[1] ± p.position_err[1]) ppm
                    δ(Y): $(p.position[2] ± p.position_err[2]) ppm
                    lw(X): $(p.lwX ± p.lwX_err) Hz
                    lw(Y): $(p.lwY ± p.lwY_err) Hz
                    """
                end
            end
        end
    end
end


function preparefigure!(state)
    mainpanel = state.fig[1:2,1]
    state.fig[1:2, 1] = GridLayout()
    # mainlayout = fig[1, 1] = GridLayout()
    # controlpanel = fig[1, 2]
    # rightcolumn = state.fig[1:2, 2] = GridLayout()
    controllayout = state.fig[1, 2] = GridLayout()
    fitplotpanel = state.fig[2, 2]
    state.fig[2, 2] = GridLayout()
    # rowsize!(rightcolumn, 1, Auto(2))
    
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

    connect!(state.setup.peakradiusX, state.controls.□Xradius.value)
    connect!(state.setup.peakradiusY, state.controls.□Yradius.value)
end


function preparefitplot!(state, panel)
    @show t = state.specdata.T
    data = lift(state.highlightedpeakindex, state.peaks) do i, P
        if i==0 || i > length(P)
            y = 0t
        else
            y = P[i].amp #.± P[i].amp_err
        end
        [Point2f(a,b) for (a,b)∈zip(t, y)]
    end
    err = lift(state.highlightedpeakindex, state.peaks) do i, P
        if i==0 || i > length(P)
            collect(0.0 * t)
        else
            # @show P[i]
            # @show P[i].amp_err
            P[i].amp_err
        end
    end

    tmin = minimum(t)
    trange = maximum(t) - tmin
    tmin -= 0.05trange
    trange *= 1.1

    limits = lift(state.highlightedpeakindex, state.peaks) do i, P
        if i==0 || i > length(P)
            ymin = 0.
            yrange = 1. 
        else
            y = P[i].amp
            ymin = minimum(y)
            yrange = maximum(y) - ymin
            ymin -= 0.05yrange
            yrange *= 1.1
        end
        (tmin, tmin+trange, ymin, ymin+yrange)
    end

    fitdata = lift(state.highlightedpeakindex, state.peaks) do i, P
        tval = LinRange(0.0, maximum(t), 50)
        if i==0 || i > length(P) || P[i].touched
            y = 0tval
        else
            y = P[i].A .* exp.(-P[i].k .* tval)
        end
        [Point2f(a,b) for (a,b)∈zip(tval, y)]
    end
    
    ax = Axis(panel, xgridvisible=false, ygridvisible=false, limits=limits, aspect=1.0)
    fitplot = lines!(ax, fitdata, color=colorant"orangered")
    errorbars!(ax, data, err, whiskerwidth=10)
    seriesplot = scatter!(ax, data)
    state.fitplots = FitPanel(ax, seriesplot, fitplot)
end


function preparecontourpanel!(state, panel)
    Xlabel = label(state.specdata.spectrum, F1Dim)
    Ylabel = label(state.specdata.spectrum, F2Dim)
    spectitle = label(state.specdata.spectrum)#[:title]
    if length(spectitle) > 60
        spectitle = spectitle[1:60] * "…"
    end

    currentZ = lift(state.controls.□plane.value) do slice
        state.specdata.Z[:,:,slice]
    end
    currentZsim = lift(state.specdata.Zfit, state.controls.□plane.value) do Zdata, slice
        Zdata[:,:,slice]
    end

    contourax = Axis(panel,
        xreversed=true, yreversed=true,
        xlabel=Xlabel, ylabel=Ylabel, title=spectitle,
        xgridvisible=false, ygridvisible=false)

    maskplot = heatmap!(contourax, state.specdata.dX, state.specdata.dY, getmask(state),
        colormap=[colorant"white", colorant"navajowhite"], inspectable=false)
    contourplot = contour!(contourax,
        state.specdata.dX,
        state.specdata.dY,
        currentZ,
        levels=state.specdata.clev,
        colormap=[colorant"red", colorant"black"],
        colorrange=(-0.001*state.specdata.σ, 0.001*state.specdata.σ),
        lowclip=:red,
        highclip=:black,
        inspectable=false)
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
    @show getpeakpositions(state)
    scatterplot = scatter!(contourax, getpeakpositions(state),
        color=getpeakcolors(state),
        markersize=15,
        marker=:circle,
        inspectable=true)
    labelsplot = text!(contourax, getpeaktext(state), offset=(10,10))

    # inspector = DataInspector(contourax)

    state.plots = PlotPanel(contourax, contourplot, scatterplot, labelsplot, maskplot, fitplot)
end

function preparecontrolpanel!(state, layout)
    layout[1,1] = state.controls.□plane#.layout
    state.controls.□plane.range[] = 1:length(state.specdata.T)
    layout[2,1] = grid!([state.controls.□contourup ;; state.controls.□contourdown])
    layout[3,1] = state.controls.□Xradius#.layout
    layout[4,1] = state.controls.□Yradius#.layout
    layout[5,1] = grid!([Label(state.fig, "peak fitting") ;; state.controls.□fitting])
    layout[6,1] = grid!([Label(state.fig, "exponential fitting") ;; state.controls.□expfitting])
    layout[7,1] = grid!([state.controls.□savelist ;; state.controls.□outputfilename])
    layout[8,1] = state.controls.□info
    
    layout.tellheight[] = true
end


function peakfit2dseries(inputfilename, taufilename, peaklistfilename=nothing)
    # load input spectrum
    specdata = loadspectrum(inputfilename, taufilename)

    # parse input peak list
    peaks = loadpeaklist(peaklistfilename, specdata)

    # setup parameters
    setup = Setup()

    # create state
    fig = Figure()
    state = State(fig, setup, specdata, peaks)

    # create figure
    preparefigure!(state)

    # set up callbacks
    # - contour levels
    on(state.controls.□contourup.clicks) do clicks
        state.specdata.clev0[] *= 1.3
    end
    on(state.controls.□contourdown.clicks) do clicks
        state.specdata.clev0[] /= 1.3
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
    # - fitting
    onany(state.peaks, state.controls.□fitting.active) do P, fitting
        if !fitting
            state.controls.□expfitting.active[] = false
        end
        fitting || return Consume(false)
        state.dragging[] && return Consume(false)

        fit!(state)
        return Consume()
    end
    # - exp fitting
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



function process_mousebutton(state, event)
    # consume input if in process of renaming a peak
    state.renaming[] && return Consume()

    if event.button == Mouse.left && event.action == Mouse.press
        plt, i = pick(state.fig.scene, Makie.mouseposition_px(state.fig.scene), 50)
        state.dragging[] = (plt == state.plots.scatterplot)
        if state.dragging[]
            @debug "event: initiate drag" i
            state.dragpeakindex = i
        end
        return Consume(state.dragging[])
    elseif event.button == Mouse.left && event.action == Mouse.release
        if state.dragging[]
            @debug "event: stop dragging"
            state.dragging[] = false
            state.clusters = updatepeakandrecluster!(state.peaks, state.dragpeakindex, state.clusters)
            @debug "event: stop dragging - done"
            notify(state.peaks)
            @debug "event: drag peak - done notifying"
            
        end
        return Consume(false)
    end
    return Consume(false)
end



function process_mouseposition(state, mousepos)
    state.renaming[] && return Consume(true)

    if state.dragging[]
        newpos = mouseposition(state.plots.axis.scene)
        state.peaks[][state.dragpeakindex].position = newpos
        state.peaks[][state.dragpeakindex].touched = true

        notify(state.peaks)
        return Consume()
    else
        plt, i = pick(state.fig.scene, Makie.mouseposition_px(state.fig.scene), 50)
        if(plt == state.plots.scatterplot)
            state.highlightedpeakindex[] = i
        end
    end
    return Consume(false)
end


    
function process_keyboardbutton(state, event)
    if state.renaming[]
        if event.action == Keyboard.press && event.key == Keyboard.enter
            state.renaming[] = false
            return Consume()
        elseif event.action == Keyboard.press && event.key == Keyboard.backspace
            if length(state.peaks[][state.activepeakindex].label) > 0
                state.peaks[][state.activepeakindex].label = state.peaks[][state.activepeakindex].label[1:end-1]
                notify(state.peaks)
                return Consume()
            end
        elseif event.action == Keyboard.press && event.key == Keyboard.escape
            # restore previous label
            state.peaks[][state.activepeakindex].label = state.oldlabel
            notify(state.peaks)

            state.renaming[] = false
            return(Consume())
        end
    else
        # not renaming - parse for peak commands
        event.action == Keyboard.press || return Consume(false)

        if event.key == Keyboard.a
            # Add marker
            newpos = mouseposition(state.plots.axis.scene)
            ix = findnearest(state.specdata.dX, newpos[1])
            iy = findnearest(state.specdata.dY, newpos[2])
            amp = state.specdata.Z[ix, iy, :]
            newpeak = PeakPseudo2D(newpos, amp, 30, 30, "X$(length(state.peaks[])+1)")

            state.clusters = addpeakandrecluster!(state.peaks, newpeak)
            notify(state.peaks)
            return Consume()
        end

        plt, i = pick(state.fig.scene, Makie.mouseposition_px(state.fig.scene), 50)
        plt == state.plots.scatterplot || return Consume(false)

        if event.key == Keyboard.d
            # Delete marker
            state.clusters = deletepeakandrecluster!(state.peaks, i, state.clusters)
            notify(state.peaks)
            return Consume()
        elseif event.key == Keyboard.r
            # Rename marker
            state.activepeakindex = i
            state.oldlabel = state.peaks[][state.activepeakindex].label # store previous label
            state.peaks[][state.activepeakindex].label = "" # clear label
            state.renaming[] = true
            state.startinglabelinput = true

            notify(state.peaks)
            return Consume()
        end
    end

    return Consume(false)
end

    
function process_unicode_input(state, character)
    if state.renaming[]
        if state.startinglabelinput && character == 'r'
            # discard 'r' character hanging over from initial keypress
            state.startinglabelinput = false
            return Consume()
        end
        state.peaks[][state.activepeakindex].label = state.peaks[][state.activepeakindex].label * character
        notify(state.peaks)
        return Consume()
    end
    return Consume(false)
end

    



function addpeakandrecluster!(peaks, newpeak)
    push!(peaks.val, newpeak)
    clusters = makeclusters(peaks)
    touchpeak!(peaks, length(peaks[]), clusters)
    
    return clusters
end

function deletepeakandrecluster!(peaks, peakindex, oldclusters)
    touchpeak!(peaks, peakindex, oldclusters)
    deleteat!(peaks.val, peakindex)
    clusters = makeclusters(peaks)
    
    return clusters
end

function updatepeakandrecluster!(peaks, peakindex, oldclusters)
    touchpeak!(peaks, peakindex, oldclusters)
    clusters = makeclusters(peaks)
    touchpeak!(peaks, peakindex, clusters)

    return clusters
end

function touchpeak!(peaks, peakindex, clusters)
    for cluster ∈ clusters
        if peakindex ∈ cluster
            for i ∈ cluster
                peaks[][i].touched = true
            end
            break
        end
    end
end



function fit!(state)
    length(state.peaks[]) > 0 || return

    expfitting = state.controls.□expfitting.active[]

    #(peaks, Zfit, clusters, dX, dY, Z, spec)
    # are there any touched peaks to be fitted?
    anytouched = any([peak.touched for peak ∈ state.peaks[]])
    if anytouched
        for cluster ∈ state.clusters
            # skip over untouched clusters
            any([state.peaks[][peakindex].touched for peakindex ∈ cluster]) || continue
            
            # fitchanged = fitcluster!(peaks, Zfit, cluster, dX, dY, Z, spec)
            fitchanged = fitcluster!(state, cluster)
            
            fitchanged && continue
            
            # if there is no change, untouch fitted peaks
            # and fit exponentials if required
            for peakindex ∈ cluster
                state.peaks[][peakindex].touched = false
                if expfitting
                    fitexp!(state, peakindex)
                end
            end
        end
        @debug "fit! done fitting"
    end

    # recalculate fitted spectrum
    state.specdata.Zfit[] = simpeaks(state) #(peaks, dX, dY, spec)
    @debug "fit! done notifying Zfit"

    if anytouched
        # update peak data
        notify(state.peaks)
        @debug "fit! done notifying peaks"
    end
end


function fitcluster!(state, cluster)
    #(peaks, Zfit, cluster, dX, dY, Z, spec)
    @debug "fitting cluster" cluster
    fitchanged = false

    # get state.setup.initial peak parameters
    npeaks = length(cluster)
    peakX = [state.peaks[][peakindex].position[1] for peakindex ∈ cluster]
    peakY = [state.peaks[][peakindex].position[2] for peakindex ∈ cluster]
    amp = [state.peaks[][peakindex].amp for peakindex ∈ cluster]
    lwX = [state.peaks[][peakindex].lwX for peakindex ∈ cluster]
    lwY = [state.peaks[][peakindex].lwY for peakindex ∈ cluster]

    # get spectrum region covering selected peaks
    Xmin = minimum(peakX) - state.setup.peakradiusX[]
    Xmax = maximum(peakX) + state.setup.peakradiusX[]
    Ymin = minimum(peakY) - state.setup.peakradiusY[]
    Ymax = maximum(peakY) + state.setup.peakradiusY[]
    @debug "ROI limits" xmin xmax ymin ymax

    # extract ROI
    iX = Xmin .≤ state.specdata.dX .≤ Xmax
    iY = Ymin .≤ state.specdata.dY .≤ Ymax
    Xroi = state.specdata.dX[iX]
    Yroi = state.specdata.dY[iY]
    T = state.specdata.T
    Zroi = state.specdata.Z[iX, iY, :]
    nT = length(T)
    
    # form mask for fitted ROIs
    mask = zeros(Bool, length(Xroi), length(Yroi))
    
    # extract masked data
    for i = 1:npeaks
        maskellipse!(mask, Xroi, Yroi, peakX[i], peakY[i], state.setup.peakradiusX[], state.setup.peakradiusY[])
    end
    maskedZ = [Zroi[:,:,i][mask] for i=1:nT]
    
    # simulation
    # allocate space for simulated spectra
    Zsim = zeros(length(Xroi), length(Yroi), nT)
    function sim(pars, getresid=false)
        # unpack parameter list
        pX = pars[:,1]
        pY = pars[:,2]
        plwX = pars[:,3]
        plwY = pars[:,4]
        
        for i=1:nT
            pamp = pars[:,4+i]
        
            # simulate spectra for each residue
            Zsim[:,:,i] = simpeaks(pamp, pX, pY, plwX, plwY, Xroi, Yroi, state.specdata.spectrum)
        end
        if getresid
            r = vcat([vec(maskedZ[i] - Zsim[:,:,i][mask]) for i=1:nT]...)
            return r
        else
            return Zsim
        end
    end

    # min/max bounds
    # peakpars = [amp ;; peakX ;; peakY ;; lwX ;; lwY]
    peakpars = zeros(npeaks, 4+nT)
    peakpars[:,1] .= peakX
    peakpars[:,2] .= peakY
    peakpars[:,3] .= lwX
    peakpars[:,4] .= lwY
    for i=1:nT
        peakpars[:,4+i] .= [a[i] for a ∈ amp]
    end
    minpars = 0peakpars
    maxpars = 0peakpars
    minpars[:,1] = peakpars[:,1] .- state.setup.maxpeakmovementX
    minpars[:,2] = peakpars[:,2] .- state.setup.maxpeakmovementY
    minpars[:,3:4] .= state.setup.minR2
    minpars[:,5:end] .= -1e9
    maxpars[:,1] = peakpars[:,1] .+ state.setup.maxpeakmovementX
    maxpars[:,2] = peakpars[:,2] .+ state.setup.maxpeakmovementY
    maxpars[:,3:4] .= state.setup.maxR2
    maxpars[:,5:end] .= 1e9

    # pack/unpack functions - fit everything, including chemical shifts
    pack(peakpars) = vec(peakpars)
    unpack(p) = reshape(p, :, 4+nT)
    
    # pack initial parameters into vector
    p0 = pack(peakpars)
    pmin = pack(minpars)
    pmax = pack(maxpars)
    p0[p0 .< pmin] .= pmin[p0 .< pmin]
    p0[p0 .> pmax] .= pmax[p0 .> pmax]

    # construct residuals function
    resid(p) = sim(unpack(p), true)

    @debug "Running fit..."
    fit = LsqFit.lmfit(resid, p0, Float64[], show_trace=false,
            maxIter=50, autodiff=:finiteforward,
            x_tol=1e-3, g_tol=1e-6,
            lower=pmin, upper=pmax)
    pfit = unpack(fit.param)
    @debug pfit

    # update peak parameters
    for i=1:npeaks
        peakindex = cluster[i]
        state.peaks[][peakindex].position = Point2f(pfit[i, 1], pfit[i, 2])
        state.peaks[][peakindex].lwX = pfit[i, 3]
        state.peaks[][peakindex].lwY = pfit[i, 4]
        state.peaks[][peakindex].amp = pfit[i, 5:end]
    end

    # try getting uncertainty
    try
        se = unpack(stderror(fit))
        for i=1:npeaks
            peakindex = cluster[i]
            state.peaks[][peakindex].position_err = Point2f(se[i, 1], se[i, 2])
            state.peaks[][peakindex].lwX_err = se[i, 3]
            state.peaks[][peakindex].lwY_err = se[i, 4]
            state.peaks[][peakindex].amp_err = se[i, 5:end]
        end
    catch
        for i=1:npeaks
            peakindex = cluster[i]
            state.peaks[][peakindex].amp_err = 0state.peaks[][peakindex].amp
            state.peaks[][peakindex].amp_err .= Inf
            state.peaks[][peakindex].position_err = Point2f(Inf, Inf)
            state.peaks[][peakindex].lwX_err = Inf
            state.peaks[][peakindex].lwY_err = Inf
        end
    end

    relchange = (pfit - peakpars) ./ peakpars
    fitchanged = any(relchange .> 0.01)
    # fitchanged = false

    
    # append cluster fit to Zfit
    # state.specdata.Zfit[][iX, iY] .+= sim(pfit)

    return fitchanged
end



function fitexp!(state)
    @debug "fitexp!"
    state.controls.□expfitting.active[] || return
    for i=1:length(state.peaks[])
        fitexp!(state, i)
    end
    @debug "fitexp! done fitting"
    # notify(state.peaks)
    # @info "fitexp! done notifying"
end



function fitexp!(state, peakindex)
    @debug "fitexp! (single peak)" peakindex
    state.controls.□expfitting.active[] || return
    state.peaks[][peakindex].touched && return

    t = state.specdata.T
    y = state.peaks[][peakindex].amp
    ye = state.peaks[][peakindex].amp_err
    A0 = maximum(abs.(y))
    idx = findfirst(abs.(y) .<= A0/2)
    k0 = 1.0 / t[idx]
    if minimum(y) < -maximum(y)
        A0 = -A0
    end

    @. model(t, p) = p[1]*exp(-p[2]*t)
    p0 = [A0, k0]
    fit = curve_fit(model, t, y, p0)
    state.peaks[][peakindex].A = fit.param[1]
    state.peaks[][peakindex].k = fit.param[2]
    try
        fiterr = stderror(fit)
        state.peaks[][peakindex].A_err = fiterr[1]
        state.peaks[][peakindex].k_err = fiterr[2]
    catch e
        state.peaks[][peakindex].A_err = Inf
        state.peaks[][peakindex].k_err = Inf
    end 
    @debug "fitexp! (single peak) done" peakindex
end



function savepeaklist(state)
    outputfilename = state.controls.□outputfilename.stored_string[]
    @info "Saving peak list to $outputfilename"

    nt = length(state.specdata.T)

    peakdata = map(state.peaks[]) do p
        if p.touched
            [p.label, p.position[1], p.position[2]]
        else
            if state.controls.□expfitting.active[]
                vcat([p.label, p.position[1], p.position[2], p.k, p.lwX, p.lwY, p.A, p.position_err[1], p.position_err[2], p.k_err, p.lwX_err, p.lwY_err, p.A_err],
                    p.amp, p.amp_err)
            else
                vcat([p.label, p.position[1], p.position[2], missing, p.lwX, p.lwY, missing, p.position_err[1], p.position_err[2], missing, p.lwX_err, p.lwY_err, missing],
                    p.amp, p.amp_err)
            end
        end
    end
    open(outputfilename; write=true) do f
        write(f, "# label   dX(ppm)   dY(ppm)   rate(s-1)   lwX(Hz)   lwY(Hz)   expamp   dX_err(ppm)   dY_err(ppm)   rate_err   lwX_err(Hz)   lwY_err(Hz)   expamp_err   ampltidues...   ampltidue_errs\n")
        writedlm(f, peakdata)
    end
end