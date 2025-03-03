"""
    TRACTExperiment <: Experiment1D

Experiment for TRACT (TROSY and anti-TROSY Correlation) analysis.

# Fields
- `specdata_trosy`: Spectral data for TROSY experiment
- `specdata_antitrosy`: Spectral data for anti-TROSY experiment
- `regions`: Observable list of integration regions
- `isfitting`: Observable indicating if fitting is active
- `τ_trosy`: Relaxation times for TROSY experiment
- `τ_antitrosy`: Relaxation times for anti-TROSY experiment
- `field_strength`: Magnetic field strength in T
- `temperature`: Sample temperature in K
- `colors`: Colors for regions
- `state`: Observable dictionary of state variables
"""
struct TRACTExperiment <: Experiment1D
    specdata_trosy::SpecData1D
    specdata_antitrosy::SpecData1D
    regions::Observable{Vector{Region}}
    isfitting::Observable{Bool}
    
    τ_trosy::Vector{Float64}
    τ_antitrosy::Vector{Float64}
    field_strength::Float64
    temperature::Float64
    
    colors::Vector{Symbol}
    state::Observable{Dict{Symbol,Observable}}
    
    current_expt::Observable{Symbol}  # :trosy or :antitrosy
    
    function TRACTExperiment(specdata_trosy, specdata_antitrosy, 
                            τ_trosy, τ_antitrosy, field_strength, temperature)
        colors = [:red, :blue, :green, :purple, :orange, :cyan, :magenta, :brown, :pink]
        
        expt = new(
            specdata_trosy,
            specdata_antitrosy,
            Observable(Region[]),
            Observable(true),
            τ_trosy,
            τ_antitrosy,
            field_strength,
            temperature,
            colors,
            Observable{Dict{Symbol,Observable}}(),
            Observable{Symbol}(:trosy)
        )
        
        expt.state[] = preparestate(expt)
        setupexptobservables!(expt)
        
        return expt
    end
end

# Define a property to access the current specdata based on the current_expt flag
function Base.getproperty(expt::TRACTExperiment, s::Symbol)
    if s === :specdata
        current = getfield(expt, :current_expt)[]
        if current === :trosy
            return getfield(expt, :specdata_trosy)
        else
            return getfield(expt, :specdata_antitrosy)
        end
    else
        return getfield(expt, s)
    end
end

"""
    tract1d()

Interactive command-line interface to start TRACT analysis.
"""
function tract1d()
    println("Current directory: $(pwd())")
    println()

    print("Enter path to TROSY experiment (i.e. Bruker experiment folder): ")
    trosy_path = readline()
    ispath(trosy_path) || throw(SystemError("No such file or directory"))

    print("Enter path to anti-TROSY experiment: ")
    antitrosy_path = readline()
    ispath(antitrosy_path) || throw(SystemError("No such file or directory"))

    return tract1d(trosy_path, antitrosy_path)
end

"""
    tract1d(trosy::String, antitrosy::String)

Start TRACT analysis on the specified experiments.
"""
function tract1d(trosy::String, antitrosy::String)
    return tract1d(loadnmr(trosy), loadnmr(antitrosy))
end

"""
    tract1d(trosy::NMRData{T,2}, antitrosy::NMRData{T,2}) where {T}

Start TRACT analysis on the given NMR data.
"""
function tract1d(trosy::NMRData{T,2}, antitrosy::NMRData{T,2}) where {T}
    # Make copies of the data
    trosy = deepcopy(trosy)
    antitrosy = deepcopy(antitrosy)
    
    # Label the datasets
    label!(trosy, "TROSY")
    label!(antitrosy, "Anti-TROSY")
    
    # Get relaxation times
    trosy = setrelaxtimes(trosy, acqus(trosy, :vdlist), "s")
    antitrosy = setrelaxtimes(antitrosy, acqus(antitrosy, :vdlist), "s")
    
    τ_trosy = data(trosy, 2)
    τ_antitrosy = data(antitrosy, 2)
    
    # Get field strength
    B0 = 2π * 1e6 * acqus(trosy, :bf1) / gyromagneticratio(H1)
    
    # Get temperature
    temp = acqus(trosy, :te)
    print("Temperature = $temp K. Press enter to confirm or type correct value (in K): ")
    response = readline()
    if length(response) > 0
        temp = tryparse(Float64, response)
    end
    
    # Prepare spectral data
    specdata_trosy = preparespecdata(trosy, TRACTExperiment)
    specdata_antitrosy = preparespecdata(antitrosy, TRACTExperiment)
    
    # Create experiment object
    expt = TRACTExperiment(
        specdata_trosy,
        specdata_antitrosy,
        τ_trosy,
        τ_antitrosy,
        B0,
        temp
    )
    
    # Launch GUI
    gui!(expt)
    
    return expt
end

"""
    preparespecdata(spec::NMRData{T,2}, ::Type{TRACTExperiment}) where T

Prepare spectral data for TRACT experiment.
"""
function preparespecdata(spec::NMRData{T,2}, ::Type{TRACTExperiment}) where T
    # Extract data
    x = [data(spec, F1Dim) for _ in 1:size(spec, 2)]
    y = [data(spec, i, :) for i in 1:size(spec, 2)]
    σ = [spec[:noise] for _ in 1:size(spec, 2)]
    
    # Create slice labels based on relaxation times
    τ_values = data(spec, 2)
    zlabels = ["τ = $(round(τ, digits=3)) s" for τ in τ_values]
    
    # Create observables for plotting
    xplot = Observable(x[1])
    yplot = Observable(y[1])
    
    return SpecData1D(spec, x, y, σ, zlabels, xplot, yplot)
end

"""
    setupexptobservables!(expt::TRACTExperiment)

Set up observables and callbacks specific to TRACT experiment.
"""
function setupexptobservables!(expt::TRACTExperiment)
    # Call the parent method first
    invoke(setupexptobservables!, Tuple{Experiment1D}, expt)
    
    # Add callback for switching between TROSY and anti-TROSY
    on(expt.current_expt) do _
        update_specdata_display!(expt)
    end
    
    return expt
end

"""
    update_specdata_display!(expt::TRACTExperiment)

Update the displayed spectrum data for TRACT experiment.
"""
function update_specdata_display!(expt::TRACTExperiment)
    # Get current slice
    slice_idx = expt.state[][:current_slice][]
    
    # Get current experiment type
    current = expt.current_expt[]
    
    # Get corresponding specdata
    specdata = current === :trosy ? expt.specdata_trosy : expt.specdata_antitrosy
    
    # Update x and y data
    if slice_idx > 0 && slice_idx <= length(specdata.x)
        specdata.xplot[] = specdata.x[slice_idx]
        specdata.yplot[] = specdata.y[slice_idx]
    end
end

"""
    setup_region_parameters!(region::Region, expt::TRACTExperiment)

Set up parameters for TRACT experiment.
"""
function setup_region_parameters!(region::Region, expt::TRACTExperiment)
    # Basic integral parameters for both TROSY and anti-TROSY
    region.parameters[:integral_trosy] = Parameter("TROSY Integral", zeros(length(expt.τ_trosy)))
    region.parameters[:integral_antitrosy] = Parameter("Anti-TROSY Integral", zeros(length(expt.τ_antitrosy)))
    region.parameters[:error_trosy] = Parameter("TROSY Error", zeros(length(expt.τ_trosy)))
    region.parameters[:error_antitrosy] = Parameter("Anti-TROSY Error", zeros(length(expt.τ_antitrosy)))
    
    # Relaxation parameters
    region.postparameters[:A_trosy] = Parameter("TROSY Amplitude", 0.0)
    region.postparameters[:R2_trosy] = Parameter("TROSY R₂", 0.0)
    region.postparameters[:A_antitrosy] = Parameter("Anti-TROSY Amplitude", 0.0)
    region.postparameters[:R2_antitrosy] = Parameter("Anti-TROSY R₂", 0.0)
    
    # TRACT parameters
    region.postparameters[:ΔR] = Parameter("ΔR", 0.0)
    region.postparameters[:ηxy] = Parameter("ηxy", 0.0)
    region.postparameters[:τc] = Parameter("τc", 0.0)
end

"""
    integrate!(region::Region, expt::TRACTExperiment)

Calculate integrals for TRACT experiment.
"""
function integrate!(region::Region, expt::TRACTExperiment)
    @debug "Integrating region $(region.label[]) for TRACT"
    
    # Integrate for TROSY
    integrate_experiment!(region, expt, :trosy)
    
    # Integrate for anti-TROSY
    integrate_experiment!(region, expt, :antitrosy)
    
    # Mark as touched
    region.touched[] = true
    
    return region.parameters[:integral_trosy].value[], region.parameters[:integral_antitrosy].value[]
end

"""
    integrate_experiment!(region::Region, expt::TRACTExperiment, expt_type::Symbol)

Helper function to integrate either TROSY or anti-TROSY data.
"""
function integrate_experiment!(region::Region, expt::TRACTExperiment, expt_type::Symbol)
    # Select the correct data
    specdata = expt_type === :trosy ? expt.specdata_trosy : expt.specdata_antitrosy
    τ_values = expt_type === :trosy ? expt.τ_trosy : expt.τ_antitrosy
    
    # Get the number of slices
    nslices = length(specdata.y)
    
    # Allocate arrays for integrals and errors
    integrals = zeros(nslices)
    errors = zeros(nslices)
    
    # Iterate through slices
    for i in 1:nslices
        # Get x and y data for this slice
        x = specdata.x[i]
        y = specdata.y[i]
        
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
    integral_key = Symbol("integral_$(expt_type)")
    error_key = Symbol("error_$(expt_type)")
    region.parameters[integral_key].value[] = integrals
    region.parameters[error_key].value[] = errors
end

"""
    fit!(region::Region, expt::TRACTExperiment)

Fit the TRACT model to the region data.
"""
function fit!(region::Region, expt::TRACTExperiment)
    @debug "Fitting TRACT for region $(region.label[])"
    
    # Get integration data
    trosy_integrals = region.parameters[:integral_trosy].value[]
    trosy_errors = region.parameters[:error_trosy].value[]
    antitrosy_integrals = region.parameters[:integral_antitrosy].value[]
    antitrosy_errors = region.parameters[:error_antitrosy].value[]
    
    # Create exponential decay model
    model_func(t, p) = p[1] * exp.(-p[2] * t)
    
    # Fit TROSY data
    p0_trosy = [maximum(trosy_integrals), 5.0]  # Initial guesses
    trosy_fit = curve_fit(model_func, expt.τ_trosy, trosy_integrals, p0_trosy, 1.0 ./ trosy_errors.^2)
    trosy_pars = coef(trosy_fit) .± stderror(trosy_fit)
    
    # Fit anti-TROSY data
    p0_antitrosy = [maximum(antitrosy_integrals), 5.0]  # Initial guesses
    antitrosy_fit = curve_fit(model_func, expt.τ_antitrosy, antitrosy_integrals, p0_antitrosy, 1.0 ./ antitrosy_errors.^2)
    antitrosy_pars = coef(antitrosy_fit) .± stderror(antitrosy_fit)
    
    # Extract parameters
    A_trosy, R2_trosy = trosy_pars
    A_antitrosy, R2_antitrosy = antitrosy_pars
    
    # Update parameters
    region.postparameters[:A_trosy].value[] = Measurements.value(A_trosy)
    region.postparameters[:A_trosy].uncertainty[] = Measurements.uncertainty(A_trosy)
    region.postparameters[:R2_trosy].value[] = Measurements.value(R2_trosy)
    region.postparameters[:R2_trosy].uncertainty[] = Measurements.uncertainty(R2_trosy)
    
    region.postparameters[:A_antitrosy].value[] = Measurements.value(A_antitrosy)
    region.postparameters[:A_antitrosy].uncertainty[] = Measurements.uncertainty(A_antitrosy)
    region.postparameters[:R2_antitrosy].value[] = Measurements.value(R2_antitrosy)
    region.postparameters[:R2_antitrosy].uncertainty[] = Measurements.uncertainty(R2_antitrosy)
    
    # Calculate TRACT parameters
    ΔR = R2_antitrosy - R2_trosy
    region.postparameters[:ΔR].value[] = Measurements.value(ΔR)
    region.postparameters[:ΔR].uncertainty[] = Measurements.uncertainty(ΔR)
    
    # Calculate ηxy (cross-correlation rate)
    ηxy = ΔR / 2
    region.postparameters[:ηxy].value[] = Measurements.value(ηxy)
    region.postparameters[:ηxy].uncertainty[] = Measurements.uncertainty(ηxy)
    
    # Mark as no longer touched and post-fitted
    region.touched[] = false
    region.postfitted[] = true
end

"""
    postfit!(region::Region, expt::TRACTExperiment)

Calculate correlation time after fitting relaxation rates.
"""
function postfit!(region::Region, expt::TRACTExperiment)
    # Calculate correlation time (τc) from ηxy
    
    # Constants
    μ0 = 4π * 1e-7  # Vacuum permeability
    γH = gyromagneticratio(H1)  # H gyromagnetic ratio
    γN = gyromagneticratio(N15)  # N gyromagnetic ratio
    ħ = 6.626e-34 / 2π  # Reduced Planck constant
    rNH = 1.02e-10  # N-H bond length
    ΔδN = 160e-6  # N CSA
    θ = 17 * π / 180  # Angle
    ωN = 2π * γN * expt.field_strength  # N Larmor frequency
    
    # Dipole-CSA cross-correlation parameters
    p = μ0 * γH * γN * ħ / (8π * sqrt(2) * rNH^3)
    c = expt.field_strength * γN * ΔδN / (3 * sqrt(2))
    f = p * c * (3*cos(θ)^2 - 1)
    
    # Get ηxy value with uncertainty
    ηxy_val = region.postparameters[:ηxy].value[][1]
    ηxy_err = region.postparameters[:ηxy].uncertainty[][1]
    ηxy = ηxy_val ± ηxy_err
    
    # Calculate τc using the analytical solution
    # This is a simplified version of the complex equation in the original code
    # More accurate than simple approximations for larger molecules
    x2 = 21952 * f^6 * ωN^6 - 3025 * f^4 * ηxy_val^2 * ωN^8 + 625 * f^2 * ηxy_val^4 * ωN^10
    x = sqrt(x2)
    y3 = 1800 * f^2 * ηxy_val * ωN^4 + 125 * ηxy_val^3 * ωN^6 + 24 * sqrt(3) * x
    y = cbrt(y3)
    
    τc = (5 * ηxy_val) / (24 * f) -
         (336 * f^2 * ωN^2 - 25 * ηxy_val^2 * ωN^4) / (24 * f * ωN^2 * y) + 
         y / (24 * f * ωN^2)
    
    # Convert to ns
    τc_ns = 1e9 * τc
    
    # Estimate uncertainty (simplified error propagation)
    τc_rel_err = ηxy_err / ηxy_val
    τc_err_ns = τc_rel_err * τc_ns
    
    # Update parameter
    region.postparameters[:τc].value[] = τc_ns
    region.postparameters[:τc].uncertainty[] = τc_err_ns
    
    region.postfitted[] = true
end

"""
    get_fit_data(region::Region, expt::TRACTExperiment)

Get data for fitting and visualization.
"""
function get_fit_data(region::Region, expt::TRACTExperiment)
    # For visualization, return both TROSY and anti-TROSY data
    # This is handled separately in get_tract_data
    return Float64[], Float64[]
end

"""
    get_tract_data(region::Region, expt::TRACTExperiment)

Get combined TROSY and anti-TROSY data for visualization.
"""
function get_tract_data(region::Region, expt::TRACTExperiment)
    isnothing(region) && return (Point2f[], [], Point2f[], [])
    
    # Get TROSY and anti-TROSY data
    trosy_x = expt.τ_trosy
    trosy_y = region.parameters[:integral_trosy].value[]
    trosy_err = region.parameters[:error_trosy].value[]
    
    antitrosy_x = expt.τ_antitrosy
    antitrosy_y = region.parameters[:integral_antitrosy].value[]
    antitrosy_err = region.parameters[:error_antitrosy].value[]
    
    # Normalize by initial intensity if fitted
    if region.postfitted[]
        I0_trosy = region.postparameters[:A_trosy].value[][1]
        I0_antitrosy = region.postparameters[:A_antitrosy].value[][1]
        
        trosy_y = trosy_y ./ I0_trosy
        trosy_err = trosy_err ./ I0_trosy
        antitrosy_y = antitrosy_y ./ I0_antitrosy
        antitrosy_err = antitrosy_err ./ I0_antitrosy
    else
        # Normalize by maximum if not fitted
        trosy_y = trosy_y ./ maximum(trosy_y)
        trosy_err = trosy_err ./ maximum(trosy_y)
        antitrosy_y = antitrosy_y ./ maximum(antitrosy_y)
        antitrosy_err = antitrosy_err ./ maximum(antitrosy_y)
    end
    
    # Create points for plot
    trosy_points = Point2f.(trosy_x, trosy_y)
    trosy_errors = [(trosy_x[i], trosy_y[i], trosy_err[i]) for i in 1:length(trosy_y)]
    
    antitrosy_points = Point2f.(antitrosy_x, antitrosy_y)
    antitrosy_errors = [(antitrosy_x[i], antitrosy_y[i], antitrosy_err[i]) for i in 1:length(antitrosy_y)]
    
    # Create fit lines if fitted
    trosy_fit_points = Point2f[]
    antitrosy_fit_points = Point2f[]
    
    if region.postfitted[]
        # Create model functions
        trosy_func(t) = exp(-region.postparameters[:R2_trosy].value[][1] * t)
        antitrosy_func(t) = exp(-region.postparameters[:R2_antitrosy].value[][1] * t)
        
        # Generate points for smooth curves
        t_max = max(maximum(trosy_x), maximum(antitrosy_x))
        t_range = range(0, t_max * 1.1, 100)
        
        trosy_fit_points = Point2f.(t_range, trosy_func.(t_range))
        antitrosy_fit_points = Point2f.(t_range, antitrosy_func.(t_range))
    end
    
    return (trosy_points, trosy_errors, trosy_fit_points, 
            antitrosy_points, antitrosy_errors, antitrosy_fit_points)
end

"""
    regioninfotext(expt::TRACTExperiment, idx)

Custom region info text for TRACT experiments.
"""
function regioninfotext(expt::TRACTExperiment, idx)
    if idx == 0
        return "No region selected"
    end
    
    region = expt.regions[][idx]
    
    if !region.postfitted[]
        return "Region: $(region.label[])\nNot fitted"
    end
    
    # Format relaxation rates
    R2_trosy = region.postparameters[:R2_trosy].value[][1]
    R2_trosy_err = region.postparameters[:R2_trosy].uncertainty[][1]
    R2_antitrosy = region.postparameters[:R2_antitrosy].value[][1]
    R2_antitrosy_err = region.postparameters[:R2_antitrosy].uncertainty[][1]
    
    # Format correlation time
    τc = region.postparameters[:τc].value[][1]
    τc_err = region.postparameters[:τc].uncertainty[][1]
    
    return "Region: $(region.label[])\n" *
           "Range: $(round(region.xstart[], digits=3)) - $(round(region.xend[], digits=3)) ppm\n" *
           "Width: $(round(width(region), digits=3)) ppm\n" *
           "\n" *
           "TROSY R₂: $(round(R2_trosy, digits=2)) ± $(round(R2_trosy_err, digits=2)) s⁻¹\n" *
           "Anti-TROSY R₂: $(round(R2_antitrosy, digits=2)) ± $(round(R2_antitrosy_err, digits=2)) s⁻¹\n" *
           "ΔR: $(round(region.postparameters[:ΔR].value[][1], digits=2)) ± $(round(region.postparameters[:ΔR].uncertainty[][1], digits=2)) s⁻¹\n" *
           "τc: $(round(τc, digits=1)) ± $(round(τc_err, digits=1)) ns"
end

"""
    experimentinfo(expt::TRACTExperiment)

Custom experiment info for TRACT experiments.
"""
function experimentinfo(expt::TRACTExperiment)
    trosy_filename = isnothing(expt.specdata_trosy.nmrdata) ? "unknown" : expt.specdata_trosy.nmrdata[:filename]
    antitrosy_filename = isnothing(expt.specdata_antitrosy.nmrdata) ? "unknown" : expt.specdata_antitrosy.nmrdata[:filename]
    
    # Format field strength
    field_T = round(expt.field_strength, digits=2)
    field_MHz = round(gyromagneticratio(H1) * expt.field_strength / (2π * 1e6), digits=1)
    
    return "Analysis type: TRACT\n" *
           "TROSY filename: $trosy_filename\n" *
           "Anti-TROSY filename: $antitrosy_filename\n" *
           "TROSY points: $(length(expt.τ_trosy))\n" *
           "Anti-TROSY points: $(length(expt.τ_antitrosy))\n" *
           "Field strength: $field_T T ($field_MHz MHz)\n" *
           "Temperature: $(round(expt.temperature)) K\n" *
           "Number of regions: $(nregions(expt))"
end