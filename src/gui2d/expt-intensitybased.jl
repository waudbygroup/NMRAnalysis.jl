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

    x::Vector{Float64}
    model::FittingModel
    visualisation::VisualisationStrategy

    function IntensityExperiment(specdata, peaks, model, xvalues=nothing, visualisation=CrossSectionVisualisation())
        if isnothing(xvalues)
            xvalues = 1.0 * collect(1:length(specdata.z))
        end
        expt = new(specdata, peaks,
                   Observable(Vector{Vector{Int}}()), # clusters
                   Observable(Vector{Bool}()), # touched
                   Observable(true), # isfitting
                   Observable(0.03; ignore_equal_values=true), # xradius
                   Observable(0.2; ignore_equal_values=true), # yradius
                   Observable{Dict}(),
                   xvalues, model, visualisation)
        setupexptobservables!(expt)
        expt.state[] = preparestate(expt)
        expt
    end
end

visualisationtype(expt::IntensityExperiment) = expt.visualisation

"""
    intensities2d(inputfilenames)

Create an intensity analysis experiment.
"""
function intensities2d(inputfilenames)
    specdata = preparespecdata(inputfilenames, IntensityExperiment)
    peaks = Observable(Vector{Peak}())
    
    expt = IntensityExperiment(
        specdata,
        peaks,
        NoFitting()
    )
    
    gui!(expt)
end


"""
    relaxation2d(inputfilenames, relaxationtimes)

Create an intensity analysis experiment.
"""
function relaxation2d(inputfilenames, relaxationtimes)
    specdata = preparespecdata(inputfilenames, IntensityExperiment)
    peaks = Observable(Vector{Peak}())

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

    # specdata.zlabels .= map(t -> "τ = $t", tau)
    
    expt = IntensityExperiment(
        specdata,
        peaks,
        ExponentialModel(),
        tau,
        ModelFitVisualisation()
    )
    
    gui!(expt)
end


"""
    recovery2d(inputfilenames, relaxationtimes)

Create an intensity analysis experiment fitted to magnetisation recovery.
"""
function recovery2d(inputfilenames, relaxationtimes)
    specdata = preparespecdata(inputfilenames, IntensityExperiment)
    peaks = Observable(Vector{Peak}())

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

    # specdata.zlabels .= map(t -> "τ = $t", tau)
    
    expt = IntensityExperiment(
        specdata,
        peaks,
        RecoveryModel(),
        tau,
        ModelFitVisualisation()
    )
    
    gui!(expt)
end



"""
    modelfit2d(inputfilenames, xvalues, equation, parameters)

Create an intensity analysis experiment fitted to a custom equation
"""
function modelfit2d(inputfilenames, xvalues, modelfunction::String, parameters::Vector{Pair{String,Float64}}, xlabel="x")
    specdata = preparespecdata(inputfilenames, IntensityExperiment)
    peaks = Observable(Vector{Peak}())

    # First handle relaxation times
    xval = Float64[]
    if xvalues isa String
        append!(xval, vec(readdlm(xvalues; comments=true)))
    elseif xvalues isa Vector
        for x in xvalues
            if x isa String
                append!(xval, vec(readdlm(x; comments=true)))
            else
                append!(xval, x)
            end
        end
    end

    # specdata.zlabels .= map(t -> "τ = $t", tau)
    model = CustomModel(modelfunction, parameters::Vector{Pair{String,Float64}}, xlabel)
    expt = IntensityExperiment(
        specdata,
        peaks,
        model,
        xval,
        ModelFitVisualisation()
    )
    
    gui!(expt)
end



# load the NMR data and prepare the SpecData object
function preparespecdata(inputfilenames, ::Type{IntensityExperiment})
    @debug "Preparing spec data for intensity experiment: $inputfilenames"

    spec, x, y, z, σ, zlabels = if inputfilenames isa String
        # load a single file
        spec, x, y, z, σ = loadspecdata(inputfilenames, IntensityExperiment)
        (SingleElementVector(spec),
            SingleElementVector(x),
            SingleElementVector(y),
            z ./ σ,
            SingleElementVector(1),
            SingleElementVector(choptitle(label(spec))))
    elseif inputfilenames isa Vector{String}
        # load multiple files
        tmp = loadspecdata.(inputfilenames, IntensityExperiment)
        spec = []
        x = []
        y = []
        z = []
        σ = []
        zlabels = []
        for t in tmp
            n = length(t[4])
            if n == 1 # z is a single slice
                push!(spec, t[1])
                push!(x, t[2])
                push!(y, t[3])
                push!(z, t[4][1])
                push!(σ, t[5])
                push!(zlabels, choptitle(label(t[1])))
            else
                append!(spec, fill(t[1], n))
                append!(x, fill(t[2], n))
                append!(y, fill(t[3], n))
                append!(z, t[4])
                append!(σ, fill(t[5], n))
                append!(zlabels, fill(choptitle(label(t[1])), n))
            end
        end
        map(MaybeVector, (spec, x, y, z ./ σ[1], σ, zlabels))
    end

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

    # Add post-parameters based on model type
    setup_post_parameters!(newpeak, expt.model)

    push!(expt.peaks[], newpeak)
    return notify(expt.peaks)
end

# No post-parameters needed for NoFitting
function setup_post_parameters!(::Peak, ::NoFitting) end

# Add post-parameters for parametric models
function setup_post_parameters!(peak::Peak, model::ParametricModel)
    for name in model.param_names
        peak.postparameters[Symbol(name)] = Parameter(name, 0.0)
    end
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

function postfit!(peak::Peak, expt::IntensityExperiment)
    postfit!(peak, expt, expt.model)
end

function get_model_data(peak, expt::IntensityExperiment)
    get_model_data(peak, expt, expt.model)
end

function peakinfotext(expt::IntensityExperiment, idx)
    if idx == 0
        return "No peak selected"
    end
    
    peak = expt.peaks[][idx]
    if !peak.postfitted[]
        return "Peak: $(peak.label[])\nNot fitted"
    end
    
    # Common peak information
    info = [
        "Peak: $(peak.label[])",
        "",
        "δX: $(peak.parameters[:x].value[][1] ± peak.parameters[:x].uncertainty[][1]) ppm",
        "δY: $(peak.parameters[:y].value[][1] ± peak.parameters[:y].uncertainty[][1]) ppm",
        "X Linewidth: $(peak.parameters[:R2x].value[][1] ± peak.parameters[:R2x].uncertainty[][1]) s⁻¹",
        "Y Linewidth: $(peak.parameters[:R2y].value[][1] ± peak.parameters[:R2y].uncertainty[][1]) s⁻¹"
    ]
    
    # Add model-specific parameters
    append!(info, model_parameter_text(peak, expt.model))
    
    join(info, "\n")
end

function experimentinfo(expt::IntensityExperiment)
    info = [
        "Analysis type: Intensity",
        "Model: $(typeof(expt.model))",
        "Filename: $(expt.specdata.nmrdata[1][:filename])",
        "Number of peaks: $(length(expt.peaks[]))",
        "Experiment title: $(expt.specdata.nmrdata[1][:title])"
    ]
    
    # Add model-specific information
    append!(info, model_info_text(expt.model, expt.x))
    
    join(info, "\n")
end

get_model_xlabel(expt::IntensityExperiment) = expt.model.xlabel
get_model_ylabel(::IntensityExperiment) = "Peak amplitude"


function slicelabel(expt::IntensityExperiment, idx)
    if length(expt.specdata.zlabels) == 1
        "Slice $idx of $(nslices(expt))"
    else
        "$(expt.specdata.zlabels[idx]) ($idx of $(nslices(expt)))"
    end
end