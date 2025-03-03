"""
    DiffusionExperiment <: Experiment1D

Experiment for analyzing diffusion data from gradient-encoded spectra.

# Fields
- `specdata`: Spectral data and metadata
- `regions`: Observable list of integration regions
- `isfitting`: Observable indicating if fitting is active
- `gradients`: Gradient strengths for each slice
- `γ`: Gyromagnetic ratio
- `δ`: Gradient pulse length
- `Δ`: Diffusion delay
- `σ`: Shape factor
- `solvent`: Solvent type (for viscosity calculation)
- `temperature`: Sample temperature in K
- `colors`: Colors for regions
- `state`: Observable dictionary of state variables
"""
struct DiffusionExperiment <: Experiment1D
    specdata::SpecData1D
    regions::Observable{Vector{Region}}
    isfitting::Observable{Bool}
    
    gradients::Vector{Float64}
    γ::Float64
    δ::Float64
    Δ::Float64
    σ::Float64
    solvent::Union{Symbol,Nothing}
    temperature::Float64
    
    colors::Vector{Symbol}
    state::Observable{Dict{Symbol,Observable}}
    
    function DiffusionExperiment(specdata, regions, gradients, γ, δ, Δ, σ, solvent, temperature)
        colors = [:red, :blue, :green, :purple, :orange, :cyan, :magenta, :brown, :pink]
        
        expt = new(
            specdata,
            regions,
            Observable(true),
            gradients,
            γ,
            δ,
            Δ,
            σ,
            solvent,
            temperature,
            colors,
            Observable{Dict{Symbol,Observable}}()
        )
        
        expt.state[] = preparestate(expt)
        setupexptobservables!(expt)
        
        return expt
    end
end

"""
    diffusion1d()

Interactive command-line interface to start diffusion analysis.
"""
function diffusion1d()
    println("Current directory: $(pwd())")
    println()

    print("Enter path to diffusion experiment (i.e. Bruker experiment folder): ")
    experiment = readline()
    ispath(experiment) || throw(SystemError("No such file or directory"))

    return diffusion1d(experiment)
end

"""
    diffusion1d(experiment::String, coherence=SQ(H1))

Start diffusion analysis on the specified experiment.
"""
function diffusion1d(experiment::String, coherence=SQ(H1))
    return diffusion1d(loadnmr(experiment), coherence)
end

"""
    diffusion1d(spec::NMRData{T,2}, coherence=SQ(H1)) where {T}

Start diffusion analysis on the given NMR data.
"""
function diffusion1d(spec::NMRData{T,2}, coherence=SQ(H1)) where {T}
    spec = deepcopy(spec)  # Work on a copy of the data
    label!(spec, "Diffusion")
    
    # 1. Get experiment parameters
    println("Parsing experiment parameters...")
    γ = gyromagneticratio(coherence)
    δ = acqus(spec, :p, 30) * 2e-6  # Gradient pulse length
    Δ = acqus(spec, :d, 20)         # Diffusion delay
    
    # Try to identify shape factor
    gpnam = acqus(spec, :gpnam, 6)
    σ = 1.0
    if length(gpnam) ≥ 4
        if gpnam[1:4] == "SMSQ"
            σ = 0.9
        elseif gpnam[1:4] == "SINE"
            σ = 0.6366
        end
    end
    
    # Get temperature and solvent
    temp = acqus(spec, :te)
    solvent_str = acqus(spec, :solvent)
    solvent = nothing
    if solvent_str == "D2O"
        solvent = :d2o
    elseif solvent_str == "H2O+D2O"
        solvent = :h2o
    end
    
    # Prompt for parameter verification
    print("Gradient pulse length δ = $(1e6*δ) μs (2*p30). Press enter to confirm or type correct value (in μs): ")
    response = readline()
    if length(response) > 0
        δ = tryparse(Float64, response) * 1e-6
    end

    print("Diffusion delay Δ = $Δ s (d20). Press enter to confirm or type correct value (in s): ")
    response = readline()
    if length(response) > 0
        Δ = tryparse(Float64, response)
    end

    print("Gradient shape factor σ = $σ (gpnam6 = $gpnam). Press enter to confirm or type correct value: ")
    response = readline()
    if length(response) > 0
        σ = tryparse(Float64, response)
    end
    
    # Get gradient list
    td = spec[:, :npoints][2]
    
    Gmax = 0.55
    print("Max. gradient strength Gmax = $Gmax T m⁻¹ (typical value for Bruker systems). Press enter to confirm or type correct value (in T m⁻¹): ")
    response = readline()
    if length(response) > 0
        Gmax = tryparse(Float64, response)
    end

    print("Enter initial gradient strength (%): ")
    response = readline()
    g1 = tryparse(Float64, response) * 0.01

    print("Enter final gradient strength (%): ")
    response = readline()
    g2 = tryparse(Float64, response) * 0.01

    print("Enter gradient ramp type ('l' for linear / 'q' for quadratic / 'e' for exponential): ")
    response = readline()
    length(response) == 1 || throw(ArgumentError("Invalid input"))
    response = lowercase(response)[1]
    if response == 'l'
        g = LinRange(g1, g2, td)
    elseif response == 'q'
        g = [g1 + (g2 - g1) * (i/(td-1))^2 for i in 0:(td-1)]
    elseif response == 'e'
        g = [g1 * (g2/g1)^(i/(td-1)) for i in 0:(td-1)]
    else
        throw(ArgumentError("Invalid input"))
    end
    
    # Convert to absolute gradient strengths
    gradients = g .* Gmax
    
    # Create experiment object
    specdata = preparespecdata(spec, DiffusionExperiment)
    
    # Create experiment object
    expt = DiffusionExperiment(
        specdata,
        Observable(Region[]),
        gradients,
        γ,
        δ,
        Δ,
        σ,
        solvent,
        temp
    )
    
    # Launch GUI
    gui!(expt)
    
    return expt
end

"""
    preparespecdata(spec::NMRData{T,2}, ::Type{DiffusionExperiment}) where T

Prepare spectral data for diffusion experiment.
"""
function preparespecdata(spec::NMRData{T,2}, ::Type{DiffusionExperiment}) where T
    # Extract data
    x = [data(spec, F1Dim) for _ in 1:size(spec, 2)]
    y = [data(spec, i, :) for i in 1:size(spec, 2)]
    σ = [spec[:noise] for _ in 1:size(spec, 2)]
    
    # Create slice labels
    zlabels = ["Gradient $(i)" for i in 1:size(spec, 2)]
    
    # Create observables for plotting
    xplot = Observable(x[1])
    yplot = Observable(y[1])
    
    return SpecData1D(spec, x, y, σ, zlabels, xplot, yplot)
end

"""
    setup_region_parameters!(region::Region, expt::DiffusionExperiment)

Set up parameters for diffusion experiment.
"""
function setup_region_parameters!(region::Region, expt::DiffusionExperiment)
    # Basic integral parameter
    region.parameters[:integral] = Parameter("Integral", zeros(length(expt.gradients)))
    region.parameters[:error] = Parameter("Error", zeros(length(expt.gradients)))
    
    # Diffusion model parameters
    region.postparameters[:A] = Parameter("Amplitude", 0.0)
    region.postparameters[:D] = Parameter("Diffusion Coefficient", 0.0)
    
    # Add hydrodynamic radius if solvent is known
    if !isnothing(expt.solvent)
        region.postparameters[:rH] = Parameter("Hydrodynamic Radius", 0.0)
    end
end

"""
    integrate!(region::Region, expt::DiffusionExperiment)

Calculate integrals for diffusion experiment.
"""
function integrate!(region::Region, expt::DiffusionExperiment)
    @debug "Integrating region $(region.label[]) for diffusion"
    
    # Get the number of slices
    nslices = length(expt.specdata.y)
    
    # Allocate array for integrals and errors
    integrals = zeros(nslices)
    errors = zeros(nslices)
    
    # Iterate through slices
    for i in 1:nslices
        # Get x and y data for this slice
        x = expt.specdata.x[i]
        y = expt.specdata.y[i]
        
        # Find points within the region
        idx = findall(x -> region.xstart[] <= x <= region.xend[], x)
        
        if !isempty(idx)
            # Calculate integral using trapezoidal rule
            dx = diff(x[idx])
            ymean = (y[idx[1:end-1]] + y[idx[2:end]]) ./ 2
            integrals[i] = sum(dx .* ymean)
            
            # Estimate error from noise
            noise_idx = max(1, idx[1] - 50):min(length(x), idx[end] + 50)
            noise_idx = setdiff(noise_idx, idx)
            if !isempty(noise_idx)
                # Use standard deviation of nearby baseline as error estimate
                noise = std(y[noise_idx])
                errors[i] = noise * sqrt(length(idx))
            else
                errors[i] = 0.01 * abs(integrals[i])  # Fallback: 1% error
            end
        end
    end
    
    # Update the parameter values
    region.parameters[:integral].value[] = integrals
    region.parameters[:error].value[] = errors
    
    # Mark as touched
    region.touched[] = true
    
    return integrals
end

"""
    fit!(region::Region, expt::DiffusionExperiment)

Fit the diffusion model to the region data.
"""
function fit!(region::Region, expt::DiffusionExperiment)
    @debug "Fitting diffusion for region $(region.label[])"
    
    # Get integration data
    integrals = region.parameters[:integral].value[]
    errors = region.parameters[:error].value[]
    
    # Create model
    model = DiffusionModel(expt.γ, expt.δ, expt.Δ, expt.σ)
    
    # Get initial parameter estimates
    I0 = maximum(integrals)
    D0 = 1.0e-10  # Initial guess for D in m²/s
    p0 = [I0, D0]
    
    # Perform fit
    fit = curve_fit(model.func, expt.gradients, integrals, p0, 1.0 ./ errors.^2)
    pfit = coef(fit)
    perr = stderror(fit)
    
    # Update parameters
    region.postparameters[:A].value[] = pfit[1]
    region.postparameters[:A].uncertainty[] = perr[1]
    region.postparameters[:D].value[] = pfit[2]
    region.postparameters[:D].uncertainty[] = perr[2]
    
    # Calculate hydrodynamic radius if solvent is known
    if !isnothing(expt.solvent)
        η = viscosity(expt.solvent, expt.temperature)
        kB = 1.38e-23  # Boltzmann constant
        
        # Calculate hydrodynamic radius in Å
        D = pfit[2] ± perr[2]
        rH = kB * expt.temperature / (6 * π * η * 0.001 * D) * 1e10
        
        # Update parameter
        region.postparameters[:rH].value[] = Measurements.value(rH)
        region.postparameters[:rH].uncertainty[] = Measurements.uncertainty(rH)
    end
    
    # Mark as no longer touched and post-fitted
    region.touched[] = false
    region.postfitted[] = true
end

"""
    postfit!(region::Region, expt::DiffusionExperiment)

Additional processing after fitting for diffusion experiment.
"""
function postfit!(region::Region, expt::DiffusionExperiment)
    # Diffusion processing is already done in fit!
    # This function exists for completeness
    region.postfitted[] = true
end

"""
    get_fit_data(region::Region, expt::DiffusionExperiment)

Get data for fitting diffusion model.
"""
function get_fit_data(region::Region, expt::DiffusionExperiment)
    # For diffusion, x is gradient strength and y is integral
    return expt.gradients, region.parameters[:integral].value[]
end

"""
    regioninfotext(expt::DiffusionExperiment, idx)

Custom region info text for diffusion experiments.
"""
function regioninfotext(expt::DiffusionExperiment, idx)
    if idx == 0
        return "No region selected"
    end
    
    region = expt.regions[][idx]
    
    if !region.postfitted[]
        return "Region: $(region.label[])\nNot fitted"
    end
    
    # Format diffusion coefficient with proper units
    D = region.postparameters[:D].value[][1]
    D_err = region.postparameters[:D].uncertainty[][1]
    D_text = "D = $(D * 1e10) ± $(D_err * 1e10) × 10⁻¹⁰ m²/s"
    
    # Add hydrodynamic radius if available
    rH_text = ""
    if haskey(region.postparameters, :rH)
        rH = region.postparameters[:rH].value[][1]
        rH_err = region.postparameters[:rH].uncertainty[][1]
        rH_text = "\nrH = $(round(rH, digits=1)) ± $(round(rH_err, digits=1)) Å"
    end
    
    return "Region: $(region.label[])\n" *
           "Range: $(round(region.xstart[], digits=3)) - $(round(region.xend[], digits=3)) ppm\n" *
           "Width: $(round(width(region), digits=3)) ppm\n" *
           "\n" *
           D_text * rH_text
end

"""
    experimentinfo(expt::DiffusionExperiment)

Custom experiment info for diffusion experiments.
"""
function experimentinfo(expt::DiffusionExperiment)
    filename = isnothing(expt.specdata.nmrdata) ? "unknown" : expt.specdata.nmrdata[:filename]
    
    # Format parameters
    δ_text = "δ = $(round(expt.δ * 1e6, digits=1)) μs"
    Δ_text = "Δ = $(round(expt.Δ * 1e3, digits=1)) ms"
    σ_text = "σ = $(expt.σ)"
    
    # Add solvent and temperature if available
    solvent_text = isnothing(expt.solvent) ? "unknown" : string(expt.solvent)
    temp_text = "$(round(expt.temperature)) K"
    
    return "Analysis type: Diffusion\n" *
           "Filename: $filename\n" *
           "Gradients: $(length(expt.gradients)) levels\n" *
           "Gradient range: $(round(minimum(expt.gradients), digits=3)) - $(round(maximum(expt.gradients), digits=3)) T/m\n" *
           "$δ_text, $Δ_text, $σ_text\n" *
           "Solvent: $solvent_text\n" *
           "Temperature: $temp_text\n" *
           "Number of regions: $(nregions(expt))"
end