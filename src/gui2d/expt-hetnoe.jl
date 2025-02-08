struct HetNOEExperiment <: FixedPeakExperiment
    specdata
    peaks
    saturation

    clusters
    touched
    isfitting

    xradius
    yradius
    state

    HetNOEExperiment(specdata, peaks, saturation) = begin
        expt = new(specdata, peaks, saturation,
            Observable(Vector{Vector{Int}}()), # clusters
            Observable(Vector{Bool}()), # touched
            Observable(true), # isfitting
            Observable(0.03, ignore_equal_values=true), # xradius
            Observable(0.2, ignore_equal_values=true), # yradius
            Observable{Dict}()
            )
        setupexptobservables!(expt)
        expt.state[] = preparestate(expt)
        expt
    end
end

# create a new hetNOE experiment
function HetNOEExperiment(planefilenames, saturation::Vector{Bool})
    specdata = preparespecdata(planefilenames, saturation, HetNOEExperiment)
    peaks = Observable(Vector{Peak}())

    HetNOEExperiment(specdata, peaks, saturation)
end



# load the NMR data and prepare the SpecData object
function preparespecdata(planefilenames, saturation, ::Type{HetNOEExperiment})
    @debug "Preparing spec data for hetNOE experiment: $planefilenames"
    spectra = map(loadnmr, planefilenames)
    x = map(spec -> data(spec, F1Dim), spectra)
    y = map(spec -> data(spec, F2Dim), spectra)
    z = map(spec -> data(spec) / scale(spec), spectra)
    σ = map(spec -> spec[:noise] / scale(spec), spectra)

    zlabels = map(sat -> sat ? "Sat" : "Ref", saturation)

    SpecData(spectra, x, y,
        z ./ σ[1],
        σ ./ σ[1],
        zlabels)
end




# implementation requirements
function addpeak!(expt::HetNOEExperiment, initialposition::Point2f, label="",
                    xradius=expt.xradius[], yradius=expt.yradius[])
    expt.state[][:total_peaks][] += 1
    if label == ""
        label = "X$(expt.state[][:total_peaks][])"
    end
    @debug "Add peak $label at $initialposition"
    newpeak = Peak(initialposition, label, xradius, yradius)
    # pars: R2x, R2y, amp
    R2x0 = MaybeVector(10.)
    R2y0 = MaybeVector(10.)
    R2x = Parameter("R2x", R2x0, minvalue=1., maxvalue=100.)
    R2y = Parameter("R2y", R2y0, minvalue=1., maxvalue=100.)
    # get initial values for amplitude
    x0, y0 = initialposition
    amp0 = map(1:nslices(expt)) do i
        ix = findnearest(expt.specdata.x[i], x0)
        iy = findnearest(expt.specdata.y[i], y0)
        expt.specdata.z[i][ix, iy]
    end
    amp = Parameter("Amplitude", amp0)
    newpeak.parameters[:R2x] = R2x
    newpeak.parameters[:R2y] = R2y
    newpeak.parameters[:amp] = amp

    newpeak.postparameters[:hetnoe] = Parameter("Heteronuclear NOE", 0.7)
    newpeak.postparameters[:amp] = Parameter("Amplitude", maximum(amp0))

    push!(expt.peaks[], newpeak)
    notify(expt.peaks)
end

function simulate!(z, peak::Peak, expt::HetNOEExperiment, xbounds=nothing, ybounds=nothing)
    n = length(z)
    for i in 1:n
        # get axis references for window functions
        xaxis = dims(expt.specdata.nmrdata[i], F1Dim)
        yaxis = dims(expt.specdata.nmrdata[i], F2Dim)
        # get axis shift values
        x = isnothing(xbounds) ? expt.specdata.x[i] : expt.specdata.x[i][xbounds[i]]
        y = isnothing(ybounds) ? expt.specdata.y[i] : expt.specdata.y[i][ybounds[i]]

        x0 = peak.parameters[:x].value[][i]
        y0 = peak.parameters[:y].value[][i]
        R2x = peak.parameters[:R2x].value[][i]
        R2y = peak.parameters[:R2y].value[][i]
        amp = peak.parameters[:amp].value[][i]
        # find indices of x and y axes within peak radius of peak position
        xi = x0 .- peak.xradius[] .≤ x .≤ x0 .+ peak.xradius[]
        yi = y0 .- peak.yradius[] .≤ y .≤ y0 .+ peak.yradius[]
        xs = x[xi]
        ys = y[yi]
        # NB. scale intensities by R2x and R2y to decouple amplitude estimation from linewidth
        zx = NMRTools.NMRBase._lineshape(getω(xaxis, x0), R2x, getω(xaxis, xs), xaxis[:window], RealLineshape())
        zy = (π^2 * amp * R2x * R2y) * NMRTools.NMRBase._lineshape(getω(yaxis, y0), R2y, getω(yaxis, ys), yaxis[:window], RealLineshape())
        z[i][xi, yi] .+= zx .* zy'
    end
end



function postfit!(peak::Peak, expt::HetNOEExperiment)
    @debug "Post-fitting peak $(peak.label)" maxlog=10
    sat = expt.saturation
    A = peak.parameters[:amp].value[] .± peak.parameters[:amp].uncertainty[]

    Iref = A[sat .== false]
    Isat = A[sat .== true]
    Iref = sum(Iref) / length(Iref)
    Isat = sum(Isat) / length(Isat)
    hetNOE = Isat / Iref
    peak.postparameters[:hetnoe].value[] .= Measurements.value(hetNOE)
    peak.postparameters[:hetnoe].uncertainty[] .= Measurements.uncertainty(hetNOE)
    peak.postparameters[:amp].value[] .= Measurements.value(Iref)
    peak.postparameters[:amp].uncertainty[] .= Measurements.uncertainty(Iref)

    peak.postfitted[] = true
end

function slicelabel(expt::HetNOEExperiment, idx)
    if expt.saturation[idx]
        "Saturated ($idx of $(nslices(expt)))"
    else
        "Reference ($idx of $(nslices(expt)))"
    end
end

function peakinfotext(expt::HetNOEExperiment, idx)
    if idx == 0
        return "No peak selected"
    end
    peak = expt.peaks[][idx]
    if peak.postfitted[]
        return "Peak: $(peak.label[])\n" *
            "HetNOE: $(peak.postparameters[:hetnoe].value[][1] ± peak.postparameters[:hetnoe].uncertainty[][1])\n" *
            "\n" *
            "δX: $(peak.parameters[:x].value[][1] ± peak.parameters[:x].uncertainty[][1]) ppm\n" *
            "δY: $(peak.parameters[:y].value[][1] ± peak.parameters[:y].uncertainty[][1]) ppm\n" *
            "Amplitude: $(peak.postparameters[:amp].value[][1] ± peak.postparameters[:amp].uncertainty[][1])\n" *
            "X Linewidth: $(peak.parameters[:R2x].value[][1] ± peak.parameters[:R2x].uncertainty[][1]) s⁻¹\n" *
            "Y Linewidth: $(peak.parameters[:R2y].value[][1] ± peak.parameters[:R2y].uncertainty[][1]) s⁻¹"
    else
        return "Peak: $(peak.label[])\n" *
            "Not fitted"
    end
end



function experimentinfo(expt::HetNOEExperiment)
    "Analysis type: Heteronuclear NOE\n" *
    "Filename: $(expt.specdata.nmrdata[1][:filename])\n" *
    "Saturation: $(join(expt.saturation, ", "))\n" *
    "Number of peaks: $(length(expt.peaks[]))\n" *
    "Experiment title: $(expt.specdata.nmrdata[1][:title])\n"
end



function completestate!(state, expt::HetNOEExperiment)
    # Set up observables for the GUI
    state[:peakplot_x] = 1:nslices(expt)

    state[:peakplot_y] = lift(state[:current_peak]) do peak
        if isnothing(peak)
            [0. for i=1:nslices(expt)]
        else
            Iref = peak.postparameters[:amp].value[][1]
            [peak.parameters[:amp].value[][i] / Iref for i=1:nslices(expt)]
        end
    end
    
    state[:peakplot_xye] = lift(state[:current_peak]) do peak
        if isnothing(peak)
            [(1.0*i,0.,0.) for i=1:nslices(expt)]
        else
            Iref = peak.postparameters[:amp].value[][1]
            [(1.0*i,
              peak.parameters[:amp].value[][i] / Iref,
              peak.parameters[:amp].uncertainty[][i] / Iref)
            for i=1:nslices(expt)]
        end
    end
    
    state[:peakplot_fit_y] = lift(state[:current_peak]) do peak
        if isnothing(peak)
            [0.]
        else
            hetnoe = peak.postparameters[:hetnoe].value[][1]
            [1, hetnoe]
        end
    end
end


"""
    # peak_plot_data(peak, expt::HetNOEExperiment)

Extract plotting data for a single peak, returning observed points, error bars,
and fit line data.
"""
# function peak_plot_data(peak, expt::HetNOEExperiment)
    # # Calculate observed data points with errors
    # t = expt.saturation
    # y = peak.parameters[:amp].value[]
    # err = peak.parameters[:amp].uncertainty[]
    # obs_points = Point2f.(t, y)
    # obs_errors = [(t[i], y[i], err[i]) for i in 1:length(t)]
    
    # # Calculate fit line
    # tpred = range(0, 1.1 * maximum(t), 100)
    # A = peak.postparameters[:amp].value[][1]
    # R = peak.postparameters[:relaxationrate].value[][1]
    # ypred = A * exp.(-R * tpred)
    # fit_points = Point2f.(tpred, ypred)
    
    # return (obs_points, obs_errors, fit_points)
# end


"""
    plot_peak!(ax, peak, saturation)

Plot a single peak's data and fit onto the given axis.
"""
function plot_peak!(panel, peak, expt::HetNOEExperiment)
    x = 1:nslices(expt)

    y = if isnothing(peak)
        [0. for i=1:nslices(expt)]
    else
        Iref = peak.postparameters[:amp].value[][1]
        [peak.parameters[:amp].value[][i] / Iref for i=1:nslices(expt)]
    end
    
    xye = if isnothing(peak)
        [(1.0*i,0.,0.) for i=1:nslices(expt)]
    else
        Iref = peak.postparameters[:amp].value[][1]
        [(1.0*i,
            peak.parameters[:amp].value[][i] / Iref,
            peak.parameters[:amp].uncertainty[][i] / Iref)
        for i=1:nslices(expt)]
    end
    
    fit_y = if isnothing(peak)
        [0.]
    else
        hetnoe = peak.postparameters[:hetnoe].value[][1]
        [1, hetnoe]
    end

    ax = Axis(panel[1,1], 
        xlabel="Experiment", 
        ylabel="Relative amplitude",
        xticks=(1:nslices(expt), expt.specdata.zlabels),
        title="$(peak.label[])")
    
    hlines!(ax, [0], linewidth=0)
    hlines!(ax, fit_y, linewidth=2, color=:red)
    errorbars!(ax, xye, whiskerwidth=10)
    barplot!(ax, x, y)

end


"""
    makepeakplot!(gui, state, expt::HetNOEExperiment)

Create interactive peak plot in GUI context.
"""
function makepeakplot!(gui, state, expt::HetNOEExperiment)
    gui[:axpeakplot] = ax = Axis(gui[:panelpeakplot][1,1], 
                                xlabel="Experiment", 
                                ylabel="Relative amplitude",
                                xticks=(1:nslices(expt), expt.specdata.zlabels))
    
    hlines!(ax, [0], linewidth=0)
    hlines!(ax, state[:peakplot_fit_y], linewidth=2, color=:red)
    errorbars!(ax, state[:peakplot_xye], whiskerwidth=10)
    barplot!(ax, state[:peakplot_x], state[:peakplot_y])

    # on(state[:peakplot_xye]) do _
    #     autolimits!(ax)
    # end
end


"""
    save_peak_plots(expt::HetNOEExperiment, folder::AbstractString)

Save individual plots for each peak in the experiment to the specified folder.
"""
function save_peak_plots!(expt::HetNOEExperiment, folder::AbstractString)
    CairoMakie.activate!()
    for peak in expt.peaks[]
        fig = Figure()
        plot_peak!(fig, peak, expt)
        
        save(joinpath(folder, "peak_$(peak.label[]).pdf"), fig)
    end
    GLMakie.activate!()
end