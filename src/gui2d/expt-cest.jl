"""
    cest2d(inputfilename)

Start interactive GUI for analyzing 2D CEST (Chemical Exchange Saturation Transfer) data.

# Arguments
- `inputfilename`: NMR data file as a processed data directory containing pseudo-3D data
                   where the first plane is the reference spectrum and subsequent planes
                   are the saturation spectra

# Example:
```julia
cest2d("path/to/expno/pdata/1")
```
"""
function cest2d(inputfilename)
    expt = CESTExperiment(inputfilename)
    gui!(expt)
end

"""
    CESTExperiment <: FixedPeakExperiment

Chemical Exchange Saturation Transfer experiment with reference and saturation spectra.

# Fields
- `specdata`: Spectral data and metadata
- `peaks`: Observable list of peaks
- `frequencies`: Vector of saturation frequencies in ppm
"""
struct CESTExperiment <: FixedPeakExperiment
    specdata::Any
    peaks::Any
    frequencies::Vector{Float64}

    clusters::Any
    touched::Any
    isfitting::Any

    xradius::Any
    yradius::Any
    state::Any

    function CESTExperiment(specdata, peaks, frequencies)
        expt = new(specdata, peaks, frequencies,
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

struct CESTVisualisation <: VisualisationStrategy end
visualisationtype(::CESTExperiment) = CESTVisualisation()

"""
    CESTExperiment(inputfilename)

Create CEST experiment from a pseudo-3D input file where the first plane is the reference
and the remaining planes are saturation spectra at different frequencies.
"""
function CESTExperiment(inputfilename)
    spec = loadnmr(inputfilename)
    
    # Extract saturation frequencies from fq3list using proper NMRTools methods
    frequencies = if haskey(acqus(spec), :fq3list)
        fq_list = acqus(spec, :fq3list)
        # Get frequencies in ppm
        getppm(fq_list, dims(spec, F2Dim))
    else
        # Fallback if fq3list is not available
        collect(range(-10.0, 10.0, length=ndims(spec, 3)))
    end
    
    # Prepare specdata
    specdata = preparespecdata(inputfilename, frequencies, CESTExperiment)
    peaks = Observable(Vector{Peak}())

    return CESTExperiment(specdata, peaks, frequencies)
end

# Load the NMR data and prepare the SpecData object
function preparespecdata(inputfilename, frequencies, ::Type{CESTExperiment})
    @debug "Preparing spec data for CEST experiment: $inputfilename"
    spec = loadnmr(inputfilename)
    x = data(spec, F1Dim)
    y = data(spec, F2Dim)
    
    # Get 3D data and normalize by scale
    raw_data = data(spec) / scale(spec)
    σ = spec[:noise] / scale(spec)
    
    # Extract slices from the 3D data
    z = eachslice(raw_data; dims=3)
    
    # Create labels for each saturation frequency
    zlabels = ["Reference"]
    for freq in frequencies[2:end]  # Skip first frequency for reference
        push!(zlabels, "$(round(freq, digits=2)) ppm")
    end

    return SpecData(SingleElementVector(spec),
            SingleElementVector(x),
            SingleElementVector(y),
            z ./ σ,
            SingleElementVector(1),
            zlabels)
end

"""Add peak to experiment, setting up type-specific parameters."""
function addpeak!(expt::CESTExperiment, initialposition::Point2f, label="",
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
        expt.specdata.z[i][ix, iy]
    end
    amp = Parameter("Amplitude", amp0)
    
    newpeak.parameters[:R2x] = R2x
    newpeak.parameters[:R2y] = R2y
    newpeak.parameters[:amp] = amp

    # Add post-parameters for CEST analysis
    newpeak.postparameters[:ref_amp] = Parameter("Reference Amplitude", 0.0)
    newpeak.postparameters[:min_intensity] = Parameter("Minimum Intensity", 0.0)
    newpeak.postparameters[:min_frequency] = Parameter("Minimum Frequency", 0.0)
    
    push!(expt.peaks[], newpeak)
    notify(expt.peaks)
end

"""Simulate single peak according to experiment type."""
function simulate!(z, peak::Peak, expt::CESTExperiment, xbounds=nothing, ybounds=nothing)
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
        zx = NMRTools.NMRBase._lineshape(getω(xaxis, x0), R2x, getω(xaxis, xs),
                                         xaxis[:window], RealLineshape())
        zy = (π^2 * amp * R2x * R2y) *
             NMRTools.NMRBase._lineshape(getω(yaxis, y0), R2y, getω(yaxis, ys),
                                         yaxis[:window], RealLineshape())
        z[i][xi, yi] .+= zx .* zy'
    end
end

"""Calculate final parameters after fitting."""
function postfit!(peak::Peak, expt::CESTExperiment)
    @debug "Post-fitting peak $(peak.label)" maxlog = 10
    
    # Get amplitudes with uncertainties
    A = peak.parameters[:amp].value[] .± peak.parameters[:amp].uncertainty[]
    
    # Reference amplitude is the first plane
    ref_amp = A[1]
    
    # Calculate relative intensities (I/I0)
    rel_intensities = A ./ ref_amp
    
    # Store reference amplitude
    peak.postparameters[:ref_amp].value[] .= Measurements.value(ref_amp)
    peak.postparameters[:ref_amp].uncertainty[] .= Measurements.uncertainty(ref_amp)
    
    # Find minimum intensity and corresponding frequency
    min_idx = findmin(Measurements.value.(rel_intensities[2:end]))[2] + 1  # +1 because we skipped first element
    min_intensity = rel_intensities[min_idx]
    min_frequency = expt.frequencies[min_idx]
    
    peak.postparameters[:min_intensity].value[] .= Measurements.value(min_intensity)
    peak.postparameters[:min_intensity].uncertainty[] .= Measurements.uncertainty(min_intensity)
    peak.postparameters[:min_frequency].value[] .= min_frequency
    
    peak.postfitted[] = true
end

"""Return descriptive text for slice idx."""
function slicelabel(expt::CESTExperiment, idx)
    if idx == 1
        "Reference"
    else
        "Saturation at $(round(expt.frequencies[idx], digits=2)) ppm ($idx of $(nslices(expt)))"
    end
end

"""Return formatted text describing peak idx."""
function peakinfotext(expt::CESTExperiment, idx)
    if idx == 0
        return "No peak selected"
    end
    
    peak = expt.peaks[][idx]
    if peak.postfitted[]
        return "Peak: $(peak.label[])\n" *
               "Min Intensity: $(peak.postparameters[:min_intensity].value[][1] ± peak.postparameters[:min_intensity].uncertainty[][1])\n" *
               "Min Frequency: $(peak.postparameters[:min_frequency].value[][1]) ppm\n" *
               "\n" *
               "δX: $(peak.parameters[:x].value[][1] ± peak.parameters[:x].uncertainty[][1]) ppm\n" *
               "δY: $(peak.parameters[:y].value[][1] ± peak.parameters[:y].uncertainty[][1]) ppm\n" *
               "Reference Amplitude: $(peak.postparameters[:ref_amp].value[][1] ± peak.postparameters[:ref_amp].uncertainty[][1])\n" *
               "X Linewidth: $(peak.parameters[:R2x].value[][1] ± peak.parameters[:R2x].uncertainty[][1]) s⁻¹\n" *
               "Y Linewidth: $(peak.parameters[:R2y].value[][1] ± peak.parameters[:R2y].uncertainty[][1]) s⁻¹"
    else
        return "Peak: $(peak.label[])\n" *
               "Not fitted"
    end
end

"""Return formatted text describing experiment."""
function experimentinfo(expt::CESTExperiment)
    return "Analysis type: CEST (Chemical Exchange Saturation Transfer)\n" *
           "Filename: $(expt.specdata.nmrdata[1][:filename])\n" *
           "Number of planes: $(nslices(expt))\n" *
           "Number of peaks: $(length(expt.peaks[]))\n" *
           "Experiment title: $(expt.specdata.nmrdata[1][:title])\n"
end

## Visualisation
function get_cest_data(peak, expt::CESTExperiment)
    @debug "getting CEST data"
    isnothing(peak) && return (Float64[], Point2f[], [(0.0, 0.0, 0.0)], 0.0)

    # X-axis will be frequency values
    x = expt.frequencies[2:end]
    
    # Get amplitudes and reference amplitude
    amp = peak.parameters[:amp].value[]
    amp_err = peak.parameters[:amp].uncertainty[]
    ref_amp = 1.0  # Reference amplitude is always 1.0
    ref_err = amp[1] > 0 ? amp_err[1] / amp[1] : 0.0
    
    # Calculate relative intensities
    y = amp[2:end] ./ amp[1]
    yerr = amp_err[2:end] ./ amp[1]
    # if error > 100% set to 100%
    yerr = min.(yerr, 1.0)
    
    # Create error tuples for plotting (without propagating reference uncertainty)
    obs_points = Point2f.(x, y)
    obs_err = collect(zip(x, y, yerr))
    
    return (x, obs_points, obs_err, ref_err)
end

function completestate!(state, expt, ::CESTVisualisation)
    @debug "completing state for CEST visualisation"
    state[:peak_plot_data] = lift(peak -> get_cest_data(peak, expt), state[:current_peak])
    state[:peak_plot_x] = lift(d -> d[1], state[:peak_plot_data])
    state[:peak_plot_obs] = lift(d -> d[2], state[:peak_plot_data])
    state[:peak_plot_err] = lift(d -> d[3], state[:peak_plot_data])
    state[:peak_plot_ref_err] = lift(d -> d[4], state[:peak_plot_data])
end

function plot_peak!(panel, peak, expt, ::CESTVisualisation)
    @debug "plotting peak for CEST visualisation"

    x, obs_points, obs_err, ref_err = get_cest_data(peak, expt)
    
    ax = Axis(panel[1, 1],
              xlabel="Saturation frequency (ppm)",
              ylabel="Relative intensity (I/I₀)")
              
    # Add reference line at y=1
    hlines!(ax, [1.0]; linewidth=1, color=:gray, linestyle=:dash)
    
    # Plot uncertainty of reference as a shaded band
    if ref_err > 0
        band!(ax, x, fill(1.0 + ref_err, length(x)), fill(1.0 - ref_err, length(x)), 
              color=(:gray, 0.2))
    end
    
    # Plot data points with error bars
    errorbars!(ax, obs_err; whiskerwidth=10)
    
    # Plot data points
    scatter!(ax, obs_points)
    
    # Connect data points with a line
    lines!(ax, obs_points)
end

function makepeakplot!(gui, state, expt, ::CESTVisualisation)
    @debug "making peak plot for CEST visualisation"
    gui[:axpeakplot] = ax = Axis(gui[:panelpeakplot][1, 1];
                                xlabel="Saturation frequency (ppm)",
                                ylabel="Relative intensity (I/I₀)",
                                title="CEST Profile")
    
    # Add reference line at y=1
    hlines!(ax, [1.0]; linewidth=1, color=:gray, linestyle=:dash)
    
    # Plot uncertainty of reference as a shaded band
    band!(ax, state[:peak_plot_x][], 
          fill(1.0 + state[:peak_plot_ref_err][], length(state[:peak_plot_x][])), 
          fill(1.0 - state[:peak_plot_ref_err][], length(state[:peak_plot_x][])), 
          color=(:gray, 0.2))
    
    
    errorbars!(ax, state[:peak_plot_err]; whiskerwidth=10)
    scatter!(ax, state[:peak_plot_obs])
    lines!(ax, state[:peak_plot_obs])
end