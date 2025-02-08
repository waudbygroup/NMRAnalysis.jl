"""
    RelaxationExperiment

NMR relaxation experiment with multiple time points.

# Fields
- `specdata`: Spectral data and metadata
- `peaks`: Observable list of peaks
- `relaxationtimes`: Vector of relaxation delay times
"""
struct RelaxationExperiment <: FixedPeakExperiment
    specdata::Any
    peaks::Any
    relaxationtimes::Any

    clusters::Any
    touched::Any
    isfitting::Any

    xradius::Any
    yradius::Any
    state::Any

    function RelaxationExperiment(specdata, peaks, relaxationtimes)
        length(relaxationtimes) == length(specdata.z) ||
        throw(ArgumentError("Number of relaxation times must match number of planes"))

        expt = new(specdata, peaks, relaxationtimes,
                   Observable(Vector{Vector{Int}}()), # clusters
                   Observable(Vector{Bool}()), # touched
                   Observable(true), # isfitting
                   Observable(0.03; ignore_equal_values=true), # xradius
                   Observable(0.2; ignore_equal_values=true), # yradius
                   Observable{Dict}())
        setupexptobservables!(expt)
        expt.state[] = preparestate(expt)
        expt
    end
end

"""
    RelaxationExperiment(experimentfiles, relaxationtimes)

Create relaxation experiment from input files and delay times
(either a list of numerical values or a file path/paths).
"""
function RelaxationExperiment(experimentfiles, relaxationtimes)
    # First handle relaxation times
    tau = Float64[]
    if relaxationtimes isa String
        append!(tau, vec(readdlm(relaxationtimes; comments=true)))
    elseif relaxationtimes isa Vector
        for t in relaxationtimes
            if t isa String
                append!(tau, vec(readdlm(t; comments=true)))
            else
                append!(tau, t)
            end
        end
    end
    zlabels = map(t -> "τ = $t", tau)

    # Then handle the input data
    spec, x, y, z, σ = if experimentfiles isa String
        # load a single file
        spec, x, y, z, σ = loadspecdata(experimentfiles, RelaxationExperiment)
        (SingleElementVector(spec), SingleElementVector(x), SingleElementVector(y), z ./ σ,
         SingleElementVector(1))
    elseif experimentfiles isa Vector{String}
        # load multiple files
        tmp = loadspecdata.(experimentfiles, RelaxationExperiment)
        spec = []
        x = []
        y = []
        z = []
        σ = []
        for t in tmp
            n = length(t[4])
            if n == 1 # z is a single slice
                push!(spec, t[1])
                push!(x, t[2])
                push!(y, t[3])
                push!(z, t[4][1])
                push!(σ, t[5])
            else
                append!(spec, fill(t[1], n))
                append!(x, fill(t[2], n))
                append!(y, fill(t[3], n))
                append!(z, t[4])
                append!(σ, fill(t[5], n))
            end
        end
        map(MaybeVector, (spec, x, y, z ./ σ[1], σ))
    end

    return RelaxationExperiment(SpecData(spec, x, y, z, σ, zlabels),
                                Observable(Vector{Peak}()),
                                tau)
end

# load the NMR data and prepare the SpecData object
function loadspecdata(inputfilename, ::Type{RelaxationExperiment})
    @debug "Loading spec data for relaxation experiment: $inputfilename"
    spec = loadnmr(inputfilename)
    x = data(spec, F1Dim)
    y = data(spec, F2Dim)

    dat = data(spec) / scale(spec)
    σ = spec[:noise] / scale(spec)

    z = if ndims(spec) == 3
        eachslice(dat; dims=3)
    else
        [dat]
    end

    return spec, x, y, z, σ
end

"""Add peak to experiment, setting up type-specific parameters."""
function addpeak!(expt::RelaxationExperiment, initialposition::Point2f, label="",
                  xradius=expt.xradius[], yradius=expt.yradius[])
    expt.state[][:total_peaks][] += 1
    if label == ""
        label = "X$(expt.state[][:total_peaks][])"
    end
    @debug "Add peak $label at $initialposition"
    newpeak = Peak(initialposition, label, xradius, yradius)
    # pars: R2x, R2y, amp
    R2x0 = MaybeVector(10.0)
    R2y0 = MaybeVector(10.0)
    R2x = Parameter("R2x", R2x0; minvalue=1.0, maxvalue=100.0)
    R2y = Parameter("R2y", R2y0; minvalue=1.0, maxvalue=100.0)
    # get initial values for amplitude
    x0, y0 = initialposition
    amp0 = map(1:nslices(expt)) do i
        ix = findnearest(expt.specdata.x[i], x0)
        iy = findnearest(expt.specdata.y[i], y0)
        return expt.specdata.z[i][ix, iy]
    end
    amp = Parameter("Amplitude", amp0)
    newpeak.parameters[:R2x] = R2x
    newpeak.parameters[:R2y] = R2y
    newpeak.parameters[:amp] = amp

    newpeak.postparameters[:relaxationrate] = Parameter("Relaxation rate",
                                                        4 / maximum(expt.relaxationtimes))
    newpeak.postparameters[:amp] = Parameter("Amplitude", maximum(amp0))

    push!(expt.peaks[], newpeak)
    return notify(expt.peaks)
end

"""Simulate single peak according to experiment type."""
function simulate!(z, peak::Peak, expt::RelaxationExperiment, xbounds=nothing,
                   ybounds=nothing)
    @debug "Simulating peak $(peak.label)" maxlog = 10
    for i in 1:nslices(expt)
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
        zx = NMRTools.NMRBase._lineshape(getω(xaxis, x0), R2x, getω(xaxis, xs),
                                         xaxis[:window], RealLineshape())
        zy = (π^2 * amp * R2x * R2y) *
             NMRTools.NMRBase._lineshape(getω(yaxis, y0), R2y, getω(yaxis, ys),
                                         yaxis[:window], RealLineshape())
        z[i][xi, yi] .+= zx .* zy'
    end
end

"""Calculate final parameters after fitting."""
function postfit!(peak::Peak, expt::RelaxationExperiment)
    @debug "Post-fitting peak $(peak.label)" maxlog = 10
    t = expt.relaxationtimes
    y = peak.parameters[:amp].value[]

    A0 = maximum(y)
    peak.postparameters[:amp].initialvalue[] .= A0
    R0 = peak.postparameters[:relaxationrate].initialvalue[][1]

    model(t, p) = p[1] * exp.(-p[2] * t)
    p0 = [A0, R0]
    fit = curve_fit(model, t, y, p0)
    pfit = coef(fit)
    perr = stderror(fit)

    peak.postparameters[:amp].value[] .= pfit[1]
    peak.postparameters[:amp].uncertainty[] .= perr[1]
    peak.postparameters[:relaxationrate].value[] .= pfit[2]
    peak.postparameters[:relaxationrate].uncertainty[] .= perr[2]

    return peak.postfitted[] = true
end

"""Return descriptive text for slice idx."""
function slicelabel(expt::RelaxationExperiment, idx)
    return "τ = $(round(expt.relaxationtimes[idx],sigdigits=3)) s ($idx of $(nslices(expt)))"
end

"""Return formatted text describing peak idx."""
function peakinfotext(expt::RelaxationExperiment, idx)
    if idx == 0
        return "No peak selected"
    end
    peak = expt.peaks[][idx]
    if peak.postfitted[]
        return "Peak: $(peak.label[])\n" *
               "Relaxation rate: $(peak.postparameters[:relaxationrate].value[][1] ± peak.postparameters[:relaxationrate].uncertainty[][1]) s⁻¹\n" *
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

"""Return formatted text describing experiment."""
function experimentinfo(expt::RelaxationExperiment)
    return "Analysis type: Relaxation experiment\n" *
           "Filename: $(expt.specdata.nmrdata[1][:filename])\n" *
           "Relaxation times: $(join(expt.relaxationtimes, ", ")) s\n" *
           "Number of peaks: $(length(expt.peaks[]))\n" *
           "Experiment title: $(expt.specdata.nmrdata[1][:title])\n"
end

function completestate!(state, expt::RelaxationExperiment)
    # Set up observables for the GUI
    state[:peakplot_obs_xy] = lift(state[:current_peak]) do peak
        return isnothing(peak) ? Point2f[] : peak_plot_data(peak, expt)[1]
    end

    state[:peakplot_obs_xye] = lift(state[:current_peak]) do peak
        return isnothing(peak) ? [(0.0, 0.0, 0.0)] : peak_plot_data(peak, expt)[2]
    end

    state[:peakplot_fit_xy] = lift(state[:current_peak]) do peak
        return isnothing(peak) ? Point2f[] : peak_plot_data(peak, expt)[3]
    end
end

"""
    peak_plot_data(peak, expt::RelaxationExperiment)

Extract plotting data for a single peak, returning observed points, error bars,
and fit line data.
"""
function peak_plot_data(peak, expt::RelaxationExperiment)
    # Calculate observed data points with errors
    t = expt.relaxationtimes
    y = peak.parameters[:amp].value[]
    err = peak.parameters[:amp].uncertainty[]
    obs_points = Point2f.(t, y)
    obs_errors = [(t[i], y[i], err[i]) for i in 1:length(t)]

    # Calculate fit line
    tpred = range(0, 1.1 * maximum(t), 100)
    A = peak.postparameters[:amp].value[][1]
    R = peak.postparameters[:relaxationrate].value[][1]
    ypred = A * exp.(-R * tpred)
    fit_points = Point2f.(tpred, ypred)

    return (obs_points, obs_errors, fit_points)
end

"""
    plot_peak!(ax, peak, relaxationtimes)

Plot a single peak's data and fit onto the given axis.
"""
function plot_peak!(ax, peak, expt::RelaxationExperiment)
    obs_points, obs_errors, fit_points = peak_plot_data(peak, expt)

    hlines!(ax, [0]; linewidth=0)
    lines!(ax, fit_points; label="Fit", color=:red)
    errorbars!(ax, obs_errors; whiskerwidth=10)
    return scatter!(ax, obs_points; label="Observed")
end

"""Set up GUI visualization."""
function makepeakplot!(gui, state, expt::RelaxationExperiment)
    gui[:axpeakplot] = ax = Axis(gui[:panelpeakplot][1, 1];
                                 xlabel="Relaxation time / s",
                                 ylabel="Amplitude")

    hlines!(ax, [0]; linewidth=0)
    lines!(ax, state[:peakplot_fit_xy]; label="Fit", color=:red)
    errorbars!(ax, state[:peakplot_obs_xye]; whiskerwidth=10)
    return scatter!(ax, state[:peakplot_obs_xy]; label="Observed")
end

"""Save publication plots for all peaks."""
function save_peak_plots!(expt::RelaxationExperiment, folder::AbstractString)
    CairoMakie.activate!()
    for peak in expt.peaks[]
        fig = Figure()
        ax = Axis(fig[1, 1];
                  xlabel="Relaxation time / s",
                  ylabel="Amplitude",
                  title="$(peak.label[])")

        plot_peak!(ax, peak, expt)
        # axislegend(ax)

        save(joinpath(folder, "peak_$(peak.label[]).pdf"), fig)
    end
    return GLMakie.activate!()
end