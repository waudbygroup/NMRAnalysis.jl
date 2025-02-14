"""
    intensities2d(inputfilenames)

Start interactive GUI for analyzing 2D NMR intensity data.

# Arguments
- `inputfilenames`: NMR data in any of these formats:
    * Single experiment directory (e.g. ".../11")
    * List of multiple experiment directories
    * List of processed data directories (e.g. ".../pdata/231")

# Example:

## Single pseudo-3d experiment
```julia
intensities2d("path/to/expno")
```
## Multiple experiments from processed data
```julia
intensities2d(
    ["path/to/expno/pdata/231", "path/to/expno/pdata/232"]
)
```
"""
function intensities2d(inputfilenames)
    expt = IntensityExperiment(inputfilenames)
    gui!(expt)
end

"""
    IntensityExperiment <: FixedPeakExperiment

Generic experiment for measurement of intensity modulations across 2D spectra.

# Fields
- `specdata`: Spectral data and metadata
- `peaks`: Observable list of peaks  
"""
struct IntensityExperiment <: FixedPeakExperiment
    specdata::Any
    peaks::Any

    clusters::Any
    touched::Any
    isfitting::Any

    xradius::Any
    yradius::Any
    state::Any

    function IntensityExperiment(specdata, peaks)
        expt = new(specdata, peaks,
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
    IntensityExperiment(inputfilenames)

Create a generic experiment for measuring 2D intensity changes from a list of 2D files or pseudo-3D spectrum.
"""
function IntensityExperiment(inputfilenames)
    specdata = preparespecdata(inputfilenames, IntensityExperiment)
    peaks = Observable(Vector{Peak}())

    return IntensityExperiment(specdata, peaks)
end

# load the NMR data and prepare the SpecData object
function preparespecdata(inputfilenames, ::Type{IntensityExperiment})
    @debug "Preparing spec data for relaxation experiment: $inputfilenames"

    spec, x, y, z, σ = if inputfilenames isa String
        # load a single file
        spec, x, y, z, σ = loadspecdata(inputfilenames, IntensityExperiment)
        (SingleElementVector(spec), SingleElementVector(x), SingleElementVector(y), z ./ σ,
         SingleElementVector(1))
    elseif inputfilenames isa Vector{String}
        # load multiple files
        tmp = loadspecdata.(inputfilenames, IntensityExperiment)
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

    zlabels = choptitle.(map(label, spec))

    return SpecData(spec, x, y, z, σ, zlabels)
end

# load the NMR data and prepare the SpecData object
function loadspecdata(inputfilename, ::Type{IntensityExperiment})
    @debug "Loading spec data for intensity experiment: $inputfilename"
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
function addpeak!(expt::IntensityExperiment, initialposition::Point2f, label="",
                  xradius=expt.xradius[], yradius=expt.yradius[])
    expt.state[][:total_peaks][] += 1
    if label == ""
        label = "X$(expt.state[][:total_peaks][])"
    end
    @debug "Add peak $label at $initialposition"
    newpeak = Peak(initialposition, label, xradius, yradius)

    # pars: R2x, R2y, amp
    R2x0 = MaybeVector(30.0)
    R2y0 = MaybeVector(15.0)
    R2x = Parameter("R2x", R2x0; minvalue=1.0, maxvalue=100.0)
    R2y = Parameter("R2y", R2y0; minvalue=1.0, maxvalue=100.0)

    # get initial values for amplitude
    x0, y0 = initialposition
    ix = findnearest(expt.specdata.x[1], x0)
    iy = findnearest(expt.specdata.y[1], y0)
    amp0 = map(1:nslices(expt)) do i
        ix = findnearest(expt.specdata.x[i], x0)
        iy = findnearest(expt.specdata.y[i], y0)
        return expt.specdata.z[i][ix, iy]
    end
    amp = Parameter("Amplitude", amp0)
    
    newpeak.parameters[:R2x] = R2x
    newpeak.parameters[:R2y] = R2y
    newpeak.parameters[:amp] = amp

    push!(expt.peaks[], newpeak)
    return notify(expt.peaks)
end

"""Simulate single peak according to experiment type."""
function simulate!(z, peak::Peak, expt::IntensityExperiment, xbounds=nothing, ybounds=nothing)
    R2x0 = peak.parameters[:R2x].value[][1]
    R2y0 = peak.parameters[:R2y].value[][1]
    amp0 = peak.parameters[:amp].value[][1]

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
        zy = (π^2 * amp * R2x0 * R2y0) *
             NMRTools.NMRBase._lineshape(getω(yaxis, y0), R2y, getω(yaxis, ys),
                                         yaxis[:window], RealLineshape())
        z[i][xi, yi] .+= zx .* zy'
    end
end

"""Calculate final parameters after fitting."""
function postfit!(peak::Peak, expt::IntensityExperiment)
    peak.postfitted[] = true
end

"""Return descriptive text for slice idx."""
function slicelabel(expt::IntensityExperiment, idx)
    return "$(expt.specdata.zlabels[idx]) ($idx of $(nslices(expt)))"
end

"""Return formatted text describing peak idx."""
function peakinfotext(expt::IntensityExperiment, idx)
    if idx == 0
        return "No peak selected"
    end
    peak = expt.peaks[][idx]
    if peak.postfitted[]
        return "Peak: $(peak.label[])\n" *
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

"""Return formatted text describing experiment."""
function experimentinfo(expt::IntensityExperiment)
    return "Analysis type: PRE experiment\n" *
           "Filename: $(expt.specdata.nmrdata[1][:filename])\n" *
           "Number of peaks: $(length(expt.peaks[]))\n" *
           "Experiment title: $(expt.specdata.nmrdata[1][:title])\n"
end

function completestate!(state, expt::IntensityExperiment)
    state[:peak_plot_data] = lift(peak -> peak_plot_data(peak, expt), state[:current_peak])
    state[:peak_plot_data_xobs] = lift(d -> d[1], state[:peak_plot_data])
    state[:peak_plot_data_yobs] = lift(d -> d[2], state[:peak_plot_data])
    state[:peak_plot_data_xfit] = lift(d -> d[3], state[:peak_plot_data])
    state[:peak_plot_data_yfit] = lift(d -> d[4], state[:peak_plot_data])
end

"""
    peak_plot_data(peak, expt::IntensityExperiment)

Extract plotting data for a single peak.
"""
function peak_plot_data(peak, expt::IntensityExperiment)
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

    for i in 1:nslices(expt)
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

"""Set up GUI visualization."""
function makepeakplot!(gui, state, expt::IntensityExperiment)
    gui[:axpeakplotX] = axX = Axis(gui[:panelpeakplot][1, 1];
                                   xlabel="dX / ppm",
                                   xreversed=true,
                                   ylabel="")
    gui[:axpeakplotY] = axY = Axis(gui[:panelpeakplot][1, 2];
                                   xlabel="dY / ppm",
                                   xreversed=true,
                                   ylabel="")

    hlines!(axX, [0]; linewidth=0)
    scatterlines!(axX, state[:peak_plot_data_xobs])
    lines!(axX, state[:peak_plot_data_xfit]; color=:red)
    hlines!(axY, [0]; linewidth=0)
    scatterlines!(axY, state[:peak_plot_data_yobs])
    lines!(axY, state[:peak_plot_data_yfit]; color=:red)
end

"""Save publication plots for all peaks."""
function save_peak_plots!(expt::IntensityExperiment, folder::AbstractString)
    @debug "Saving peak plots to $folder"
    CairoMakie.activate!()

    for peak in expt.peaks[]
        xobs, yobs, xfit, yfit = peak_plot_data(peak, expt)

        fig = Figure()

        axX = Axis(fig[1, 1];
                   xlabel="δX / ppm",
                   xreversed=true,
                   ylabel="")
        axY = Axis(fig[1, 2];
                   xlabel="δY / ppm",
                   xreversed=true,
                   ylabel="")
        @debug "Axes created"
        @debug "Plotting data" xobs yobs xfit yfit

        hlines!(axX, [0]; linewidth=0)
        scatterlines!(axX, xobs)
        lines!(axX, xfit; color=:red)
        hlines!(axY, [0]; linewidth=0)
        scatterlines!(axY, yobs)
        lines!(axY, yfit; color=:red)
        @debug "Data plotted"

        save(joinpath(folder, "peak_$(peak.label[]).pdf"), fig)
    end

    GLMakie.activate!()
end