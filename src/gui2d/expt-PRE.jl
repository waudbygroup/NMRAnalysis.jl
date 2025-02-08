struct PREExperiment <: FixedPeakExperiment
    specdata
    peaks
    paramagnetic_concs
    expttype
    Trelax

    clusters
    touched
    isfitting

    xradius
    yradius
    state

    PREExperiment(specdata, peaks, paramagnetic_concs, expttype, Trelax) = begin
        expt = new(specdata, peaks, paramagnetic_concs, expttype, Trelax,
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

# create a new relaxation experiment
function PREExperiment(inputfilenames, paramagnetic_concs, exptexperimenttype, Trelax)
    exptexperimenttype in [:hsqc, :hmqc] || throw(ArgumentError("Experiment type must be :hsqc or :hmqc"))

    specdata = preparespecdata(inputfilenames, paramagnetic_concs, PREExperiment)
    peaks = Observable(Vector{Peak}())

    PREExperiment(specdata, peaks, paramagnetic_concs, exptexperimenttype, Trelax)
end



# load the NMR data and prepare the SpecData object
function preparespecdata(inputfilenames, paramagnetic_concs, ::Type{PREExperiment})
    @debug "Preparing spec data for relaxation experiment: $inputfilenames"

    spectra = map(loadnmr, inputfilenames)
    x = map(spec -> data(spec, F1Dim), spectra)
    y = map(spec -> data(spec, F2Dim), spectra)
    z = map(spec -> data(spec) / scale(spec), spectra)
    σ = map(spec -> spec[:noise] / scale(spec), spectra)

    zlabels = map(paramagnetic_concs) do conc
        if conc == 0
            "Diamagnetic"
        else
            "Paramagnetic (conc = $conc)"
        end
    end

    SpecData(spectra, x, y,
        z ./ σ[1],
        σ ./ σ[1],
        zlabels)
end




# implementation requirements
function addpeak!(expt::PREExperiment, initialposition::Point2f, label="",
                    xradius=expt.xradius[], yradius=expt.yradius[])
    expt.state[][:total_peaks][] += 1
    if label == ""
        label = "X$(expt.state[][:total_peaks][])"
    end
    @debug "Add peak $label at $initialposition"
    newpeak = Peak(initialposition, label, xradius, yradius)

    # pars: R2x, R2y, amp
    R2x0 = MaybeVector(30.)
    R2y0 = MaybeVector(15.)
    R2x = Parameter("R2x", R2x0, minvalue=1., maxvalue=100.)
    R2y = Parameter("R2y", R2y0, minvalue=1., maxvalue=100.)

    # get initial values for amplitude
    x0, y0 = initialposition
    ix = findnearest(expt.specdata.x[1], x0)
    iy = findnearest(expt.specdata.y[1], y0)
    amp0 = expt.specdata.z[1][ix, iy]
    amp0 = MaybeVector(amp0)
    amp = Parameter("Amplitude", amp0)

    # PRE
    Γ = Parameter("PRE", MaybeVector(10.), minvalue=0., maxvalue=200.)

    newpeak.parameters[:R2x] = R2x
    newpeak.parameters[:R2y] = R2y
    newpeak.parameters[:amp] = amp
    newpeak.parameters[:PRE] = Γ

    newpeak.postparameters[:PRE] = Parameter("PRE", 0.0)

    push!(expt.peaks[], newpeak)
    notify(expt.peaks)
end

function simulate!(z, peak::Peak, expt::PREExperiment)
    R2x0 = peak.parameters[:R2x].value[][1]
    R2y0 = peak.parameters[:R2y].value[][1]
    amp0 = peak.parameters[:amp].value[][1]
    PRE = peak.parameters[:PRE].value[][1]

    n = length(z)
    for i in 1:n
        # get axis references for window functions
        xaxis = dims(expt.specdata.nmrdata[i], F1Dim)
        yaxis = dims(expt.specdata.nmrdata[i], F2Dim)
        # get axis shift values
        x = data(xaxis)
        y = data(yaxis)

        x0 = peak.parameters[:x].value[][i]
        y0 = peak.parameters[:y].value[][i]

        # apply PRE to linewidths and amplitude
        R2x = R2x0 + PRE * expt.paramagnetic_concs[i]
        if expt.expttype == :hmqc
            R2y = R2y0 + PRE * expt.paramagnetic_concs[i]
        else
            R2y = R2y0
        end
        amp = amp0 * exp(-PRE * expt.paramagnetic_concs[i] * expt.Trelax)

        # find indices of x and y axes within peak radius of peak position
        xi = x0 .- peak.xradius[] .≤ x .≤ x0 .+ peak.xradius[]
        yi = y0 .- peak.yradius[] .≤ y .≤ y0 .+ peak.yradius[]
        xs = x[xi]
        ys = y[yi]
        # NB. scale intensities by R2x and R2y to decouple amplitude estimation from linewidth
        zx = NMRTools.NMRBase._lineshape(getω(xaxis, x0), R2x, getω(xaxis, xs), xaxis[:window], RealLineshape())
        zy = (π^2 * amp * R2x0 * R2y0) * NMRTools.NMRBase._lineshape(getω(yaxis, y0), R2y, getω(yaxis, ys), yaxis[:window], RealLineshape())
        z[i][xi, yi] .+= zx .* zy'
    end
end



function postfit!(peak::Peak, expt::PREExperiment)
    peak.postparameters[:PRE].uncertainty[] .= peak.parameters[:PRE].uncertainty[]
    peak.postparameters[:PRE].value[] .= peak.parameters[:PRE].value[]
    peak.postfitted[] = true
end

function slicelabel(expt::PREExperiment, idx)
    "$(expt.specdata.zlabels[idx]) ($idx of $(nslices(expt)))"
end

function peakinfotext(expt::PREExperiment, idx)
    if idx == 0
        return "No peak selected"
    end
    peak = expt.peaks[][idx]
    if peak.postfitted[]
        return "Peak: $(peak.label[])\n" *
            "PRE: $(peak.parameters[:PRE].value[][1] ± peak.parameters[:PRE].uncertainty[][1]) s⁻¹ [conc⁻¹]\n" *
            "\n" *
            "δX: $(peak.parameters[:x].value[][1] ± peak.parameters[:x].uncertainty[][1]) ppm\n" *
            "δY: $(peak.parameters[:y].value[][1] ± peak.parameters[:y].uncertainty[][1]) ppm\n" *
            "Amplitude: $(peak.parameters[:amp].value[][1] ± peak.parameters[:amp].uncertainty[][1])\n" *
            "X Linewidth: $(peak.parameters[:R2x].value[][1] ± peak.parameters[:R2x].uncertainty[][1]) s⁻¹\n" *
            "Y Linewidth: $(peak.parameters[:R2y].value[][1] ± peak.parameters[:R2y].uncertainty[][1]) s⁻¹"
    else
        return "Peak: $(peak.label[])\n" *
            "Not fitted"
    end
end



function experimentinfo(expt::PREExperiment)
    "Analysis type: PRE experiment\n" *
    "Filename: $(expt.specdata.nmrdata[1][:filename])\n" *
    "PRE agent concs: $(join(expt.paramagnetic_concs, ", "))\n" *
    "Experiment type: $(expt.expttype==:hsqc ? "HSQC" : "HMQC")\n" *
    "Relaxation time: $(expt.Trelax)\n" *
    "Number of peaks: $(length(expt.peaks[]))\n" *
    "Experiment title: $(expt.specdata.nmrdata[1][:title])\n"
end



function completestate!(state, expt::PREExperiment)
    state[:peak_plot_data] = lift(peak -> peak_plot_data(peak, expt), state[:current_peak])
    state[:peak_plot_data_xobs] = lift(d ->d[1], state[:peak_plot_data])
    state[:peak_plot_data_yobs] = lift(d -> d[2], state[:peak_plot_data])
    state[:peak_plot_data_xfit] = lift(d -> d[3], state[:peak_plot_data])
    state[:peak_plot_data_yfit] = lift(d -> d[4], state[:peak_plot_data])
end

function flatten_with_nan_separator(vectors::Vector{Vector{Point2f}})
    isempty(vectors) && return Point2f[]
    
    separator = Point2f(NaN, NaN)
    result = reduce(vectors[2:end]; init=vectors[1]) do acc, subvector
        vcat(acc, [separator], subvector)
    end
    
    # Remove the trailing separator if it exists
    return length(result) > 0 ? result[1:end-1] : result
end


"""
    peak_plot_data(peak, expt::PREExperiment)

Extract plotting data for a single peak.
"""
function peak_plot_data(peak, expt::PREExperiment)
    @debug "Preparing peak plot data"

    if isnothing(peak)
        return (Point2f[], Point2f[], Point2f[], Point2f[])
    end
    
    # we want to plot x and y cross-sections of the observed and fitted spectra
    # i.e. a list with Vector{Point2f} cross-sections for each slice - xobs, yobs, xfit, yfit
    xobs = Vector{Point2f}[]
    yobs = Vector{Point2f}[]
    xfit = Vector{Point2f}[]
    yfit = Vector{Point2f}[]

    for i=1:nslices(expt)
        x = expt.specdata.x[i]
        y = expt.specdata.y[i]

        x0 = peak.parameters[:x].value[][i]
        y0 = peak.parameters[:y].value[][i]

        # find indices of x and y axes within peak radius of peak position
        ix = x0 .- peak.xradius[] .≤ x .≤ x0 .+ peak.xradius[]
        iy = y0 .- peak.yradius[] .≤ y .≤ y0 .+ peak.yradius[]
        ix0 = findnearest(x, x0)
        iy0 = findnearest(y, y0)
        xs = x[ix]
        ys = y[iy]

        xo = Point2f.(xs, expt.specdata.z[i][ix, iy0])
        yo = Point2f.(ys, expt.specdata.z[i][ix0, iy])
        xf = Point2f.(xs, expt.specdata.zfit[][i][ix, iy0])
        yf = Point2f.(ys, expt.specdata.zfit[][i][ix0, iy])

        push!(xobs, xo)
        push!(yobs, yo)
        push!(xfit, xf)
        push!(yfit, yf)
    end
    
    @debug "Peak plot data prepared"
    return (flatten_with_nan_separator(xobs),
        flatten_with_nan_separator(yobs),
        flatten_with_nan_separator(xfit),
        flatten_with_nan_separator(yfit))
end


"""
    makepeakplot!(gui, state, expt::PREExperiment)

Create interactive peak plot in GUI context.
"""
function makepeakplot!(gui, state, expt::PREExperiment)
    gui[:axpeakplotX] = axX = Axis(gui[:panelpeakplot][1,1], 
                                xlabel="dX / ppm",
                                xreversed=true,
                                ylabel="")
    gui[:axpeakplotY] = axY = Axis(gui[:panelpeakplot][1,2], 
                                xlabel="dY / ppm", 
                                xreversed=true,
                                ylabel="")
    
    hlines!(axX, [0], linewidth=0)
    scatterlines!(axX, state[:peak_plot_data_xobs])
    lines!(axX, state[:peak_plot_data_xfit], color=:red)
    hlines!(axY, [0], linewidth=0)
    scatterlines!(axY, state[:peak_plot_data_yobs])
    lines!(axY, state[:peak_plot_data_yfit], color=:red)
end


"""
    save_peak_plots(expt::PREExperiment, folder::AbstractString)

Save individual plots for each peak in the experiment to the specified folder.
"""
function save_peak_plots!(expt::PREExperiment, folder::AbstractString)
    @debug "Saving peak plots to $folder"
    CairoMakie.activate!()

    for peak in expt.peaks[]
        xobs, yobs, xfit, yfit = peak_plot_data(peak, expt)

        fig = Figure()

        axX = Axis(fig[1,1], 
            xlabel="δX / ppm",
            xreversed=true,
            ylabel="")
        axY = Axis(fig[1,2], 
            xlabel="δY / ppm", 
            xreversed=true,
            ylabel="")
        @debug "Axes created"
        @debug "Plotting data" xobs yobs xfit yfit

        hlines!(axX, [0], linewidth=0)
        scatterlines!(axX, xobs)
        lines!(axX, xfit, color=:red)
        hlines!(axY, [0], linewidth=0)
        scatterlines!(axY, yobs)
        lines!(axY, yfit, color=:red)
        @debug "Data plotted"
        
        save(joinpath(folder, "peak_$(peak.label[]).pdf"), fig)
    end

    GLMakie.activate!()
end