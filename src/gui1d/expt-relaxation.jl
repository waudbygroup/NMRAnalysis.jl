"""
    RelaxationExperiment <: Experiment1D

Experiment for general relaxation data from time-series spectra.

# Fields
- `specdata`: Spectral data and metadata
- `regions`: Observable list of integration regions
- `isfitting`: Observable indicating if fitting is active
- `τ_values`: Relaxation times for each slice
- `model`: The fitting model to use
- `colors`: Colors for regions
- `state`: Observable dictionary of state variables
"""
struct RelaxationExperiment <: Experiment1D
    specdata::SpecData1D
    regions::Observable{Vector{Region}}
    isfitting::Observable{Bool}
    
    τ_values::Vector{Float64}
    model::ParametricModel
    
    colors::Vector{Symbol}
    state::Observable{Dict{Symbol,Observable}}
    
    function RelaxationExperiment(specdata, regions, τ_values, model)
        colors = [:red, :blue, :green, :purple, :orange, :cyan, :magenta, :brown, :pink]
        
        expt = new(
            specdata,
            regions,
            Observable(true),
            τ_values,
            model,
            colors,
            Observable{Dict{Symbol,Observable}}()
        )
        
        expt.state[] = preparestate(expt)
        setupexptobservables!(expt)
        
        return expt
    end
end

"""
    relaxation1d(inputfilenames, relaxationtimes, model=:exponential)

Create a relaxation analysis experiment.

# Arguments
- `inputfilenames`: Path to NMR data (pseudo-2D) or list of paths (multiple 1D files)
- `relaxationtimes`: List of relaxation times or path to a text file with times
- `model`: Model type (:exponential, :recovery, or :custom)

# Examples
```julia
# Using a pseudo-2D file with exponential decay
relaxation1d("path/to/t2_experiment", [0.01, 0.02, 0.03, 0.04, 0.05])

# Using multiple 1D files with recovery model
relaxation1d(["exp1.fid", "exp2.fid", "exp3.fid"], [0.1, 0.5, 1.0], :recovery)
```
"""
function relaxation1d(inputfilenames, relaxationtimes, model=:exponential)
    # Load data
    specdata = preparespecdata(inputfilenames, RelaxationExperiment)
    
    # Process relaxation times
    τ_values = Float64[]
    if relaxationtimes isa String
        append!(τ_values, vec(readdlm(relaxationtimes; comments=true)))
    elseif relaxationtimes isa Vector
        for t in relaxationtimes
            if t isa String
                append!(τ_values, vec(readdlm(t; comments=true)))
            else
                append!(τ_values, t)
            end
        end
    end
    
    # Create fitting model based on type
    fitting_model = if model == :exponential
        ExponentialModel()
    elseif model == :recovery
        RecoveryModel()
    elseif model == :custom
        # Default custom model (y = A * exp(-R * x))
        CustomModel("A * exp(-R * x)", [
            "A" => 1.0, 
            "R" => 1.0
        ], "Time / s")
    else
        error("Unsupported model type: $model")
    end
    
    # Create experiment
    expt = RelaxationExperiment(
        specdata,
        Observable(Region[]),
        τ_values,
        fitting_model
    )
    
    # Launch GUI
    gui!(expt)
    
    return expt
end

"""
    modelfit1d(inputfilenames, xvalues, modelfunction, parameters, xlabel="x")

Create a custom model fit experiment.

# Arguments
- `inputfilenames`: Path to NMR data (pseudo-2D) or list of paths (multiple 1D files)
- `xvalues`: List of x values or path to a text file with values
- `modelfunction`: String representation of the model function (e.g., "A * exp(-R * x)")
- `parameters`: Vector of parameter name/initial value pairs
- `xlabel`: Label for the x-axis

# Example
```julia
modelfit1d("path/to/data", [1, 2, 3, 4, 5],
           "A * exp(-B * x) + C",
           ["A" => 1.0, "B" => 0.5, "C" => 0.1],
           "Time / s")
```
"""
function modelfit1d(inputfilenames, xvalues, modelfunction::String, 
                   parameters::Vector{Pair{String,Float64}}, xlabel="x")
    # Load data
    specdata = preparespecdata(inputfilenames, RelaxationExperiment)
    
    # Process x values
    x_values = Float64[]
    if xvalues isa String
        append!(x_values, vec(readdlm(xvalues; comments=true)))
    elseif xvalues isa Vector
        for x in xvalues
            if x isa String
                append!(x_values, vec(readdlm(x; comments=true)))
            else
                append!(x_values, x)
            end
        end
    end
    
    # Create custom model
    model = CustomModel(modelfunction, parameters, xlabel)
    
    # Create experiment
    expt = RelaxationExperiment(
        specdata,
        Observable(Region[]),
        x_values,
        model
    )
    
    # Launch GUI
    gui!(expt)
    
    return expt
end

"""
    setup_region_parameters!(region::Region, expt::RelaxationExperiment)

Set up parameters for relaxation experiment.
"""
function setup_region_parameters!(region::Region, expt::RelaxationExperiment)
    # Basic integral parameter
    region.parameters[:integral] = Parameter("Integral", zeros(length(expt.τ_values)))
    region.parameters[:error] = Parameter("Error", zeros(length(expt.τ_values)))
    
    # Add model-specific parameters
    for name in expt.model.param_names
        region.postparameters[Symbol(name)] = Parameter(name, 0.0)
    end
end

"""
    integrate!(region::Region, expt::RelaxationExperiment)

Calculate integrals for relaxation experiment.
"""
function integrate!(region::Region, expt::RelaxationExperiment)
    @debug "Integrating region $(region.label[]) for relaxation"
    
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
    fit!(region::Region, expt::RelaxationExperiment)

Fit the relaxation model to the region data.
"""
function fit!(region::Region, expt::RelaxationExperiment)
    @debug "Fitting relaxation for region $(region.label[])"
    
    # Get integration data
    integrals = region.parameters[:integral].value[]
    errors = region.parameters[:error].value[]
    
    # Get x values (relaxation times)
    x = expt.τ_values
    
    # Get initial parameter estimates
    p0 = collect(values(estimate_parameters(x, integrals, expt.model)))
    
    # Perform fit
    fit = curve_fit(expt.model.func, x, integrals, p0, 1.0 ./ errors.^2)
    pfit = coef(fit)
    perr = stderror(fit)
    
    # Update parameters
    for (i, name) in enumerate(expt.model.param_names)
        region.postparameters[Symbol(name)].value[] = pfit[i]
        region.postparameters[Symbol(name)].uncertainty[] = perr[i]
    end
    
    # Mark as no longer touched and post-fitted
    region.touched[] = false
    region.postfitted[] = true
end

"""
    postfit!(region::Region, expt::RelaxationExperiment)

Process results after fitting for relaxation experiment.
"""
function postfit!(region::Region, expt::RelaxationExperiment)
    # For standard relaxation models, no additional post-processing is needed
    region.postfitted[] = true
end

"""
    get_fit_data(region::Region, expt::RelaxationExperiment)

Get data for fitting relaxation model.
"""
function get_fit_data(region::Region, expt::RelaxationExperiment)
    # For relaxation, x is tau values and y is integral
    return expt.τ_values, region.parameters[:integral].value[]
end

"""
    regioninfotext(expt::RelaxationExperiment, idx)

Custom region info text for relaxation experiments.
"""
function regioninfotext(expt::RelaxationExperiment, idx)
    if idx == 0
        return "No region selected"
    end
    
    region = expt.regions[][idx]
    
    if !region.postfitted[]
        return "Region: $(region.label[])\nNot fitted"
    end
    
    # Basic region info
    info = [
        "Region: $(region.label[])",
        "Range: $(round(region.xstart[], digits=3)) - $(round(region.xend[], digits=3)) ppm",
        "Width: $(round(width(region), digits=3)) ppm",
        ""
    ]
    
    # Add parameter values
    for (name, param) in region.postparameters
        value = param.value[][1]
        uncertainty = param.uncertainty[][1]
        
        # Format based on parameter type
        if string(name) == "R"
            # For rate constants, also show time constant
            rate = value
            rate_err = uncertainty
            time_const = 1.0 / rate
            time_const_err = rate_err / (rate^2)
            
            push!(info, "R: $(round(rate, digits=3)) ± $(round(rate_err, digits=3)) s⁻¹")
            push!(info, "T: $(round(time_const * 1000, digits=1)) ± $(round(time_const_err * 1000, digits=1)) ms")
        else
            push!(info, "$(param.label): $(round(value, digits=4)) ± $(round(uncertainty, digits=4))")
        end
    end
    
    return join(info, "\n")
end

"""
    experimentinfo(expt::RelaxationExperiment)

Custom experiment info for relaxation experiments.
"""
function experimentinfo(expt::RelaxationExperiment)
    filename = isnothing(expt.specdata.nmrdata) ? "unknown" : expt.specdata.nmrdata[:filename]
    
    model_type = typeof(expt.model)
    model_name = if model_type <: ExponentialModel
        "Exponential decay (A * exp(-R * t))"
    elseif model_type <: RecoveryModel
        "Recovery (A * (1 - C * exp(-R * t)))"
    elseif model_type <: CustomModel
        "Custom model"
    else
        string(model_type)
    end
    
    return "Analysis type: Relaxation\n" *
           "Model: $model_name\n" *
           "Filename: $filename\n" *
           "Time points: $(length(expt.τ_values))\n" *
           "Time range: $(minimum(expt.τ_values)) - $(maximum(expt.τ_values)) $(expt.model.xlabel)\n" *
           "Number of regions: $(nregions(expt))"
end

"""
    update_model_visualization!(g, state, expt::RelaxationExperiment, region)

Update visualization for relaxation experiment.
"""
function update_model_visualization!(g, state, expt::RelaxationExperiment, region)
    # Get model data
    obs_points, obs_errors, fit_points = get_model_data(region, expt, expt.model)
    
    # Update plot data
    g[:pltdata][] = obs_points
    g[:plterrors][] = obs_errors
    g[:pltfit][] = fit_points
    
    # Update axis labels
    g[:axplot].xlabel = expt.model.xlabel
    g[:axplot].ylabel = "Integrated signal"
end

"""
    create_fit_plot!(fig, panel, expt::RelaxationExperiment, region::Region)

Create relaxation-specific fit plot.
"""
function create_fit_plot!(fig, panel, expt::RelaxationExperiment, region::Region)
    # Create axis
    ax = Axis(panel[1, 1];
             xlabel=expt.model.xlabel,
             ylabel="Integrated signal",
             title="Relaxation Fit")
    
    if !region.postfitted[]
        text!(ax, 0.5, 0.5;
             text="No fit results available",
             align=(:center, :center),
             space=:relative)
        return
    end
    
    # Get data
    integrals = region.parameters[:integral].value[]
    errors = region.parameters[:error].value[]
    x = expt.τ_values
    
    # Create scatter plot of data
    scatter!(ax, x, integrals;
            color=region.color[],
            markersize=8)
    
    # Add error bars
    for i in 1:length(x)
        lines!(ax, [x[i], x[i]], 
              [integrals[i] - errors[i], integrals[i] + errors[i]];
              color=region.color[])
    end
    
    # Create fit line
    x_range = range(0, maximum(x) * 1.1, 100)
    parameters = [region.postparameters[Symbol(name)].value[][1] for name in expt.model.param_names]
    y_fit = expt.model.func(x_range, parameters)
    
    lines!(ax, x_range, y_fit;
          color=:red,
          linewidth=2)
    
    # Add fit results as text
    result_text = ""
    for name in expt.model.param_names
        param = region.postparameters[Symbol(name)]
        result_text *= "$(param.label): $(round(param.value[][1], digits=4)) ± $(round(param.uncertainty[][1], digits=4))\n"
    end
    
    # For rate constants, also show time constant
    if haskey(region.postparameters, :R)
        rate = region.postparameters[:R].value[][1]
        rate_err = region.postparameters[:R].uncertainty[][1]
        time_const = 1.0 / rate
        time_const_err = rate_err / (rate^2)
        
        result_text *= "T (1/R): $(round(time_const * 1000, digits=1)) ± $(round(time_const_err * 1000, digits=1)) ms\n"
    end
    
    text!(ax, 0.05, 0.95;
         text=result_text,
         align=(:left, :top),
         space=:relative)
end

"""
    save_experiment_results!(expt::RelaxationExperiment, folder)

Save additional results specific to relaxation experiments.
"""
function save_experiment_results!(expt::RelaxationExperiment, folder)
    # Save detailed relaxation parameters
    filepath = joinpath(folder, "relaxation-parameters.txt")
    
    open(filepath, "w") do f
        println(f, "# Relaxation Parameters")
        println(f, "# Generated by NMRAnalysis.jl")
        println(f, "#")
        
        # Write model information
        model_type = typeof(expt.model)
        model_name = if model_type <: ExponentialModel
            "Exponential decay (A * exp(-R * t))"
        elseif model_type <: RecoveryModel
            "Recovery (A * (1 - C * exp(-R * t)))"
        elseif model_type <: CustomModel
            "Custom model"
        else
            string(model_type)
        end
        
        println(f, "Model: $model_name")
        println(f, "#")
        
        # Write time point list
        println(f, "# Time points ($(expt.model.xlabel))")
        for (i, t) in enumerate(expt.τ_values)
            println(f, "$i\t$t")
        end
        println(f, "#")
        
        # Write integration data for each region
        println(f, "# Region integration data")
        println(f, "# Format: time_idx\ttime\tregion1\tregion1_error\tregion2\t...")
        
        # Header with region labels
        header = ["idx", "time"]
        for region in expt.regions[]
            push!(header, region.label[], "$(region.label[])_error")
        end
        println(f, "# ", join(header, "\t"))
        
        # Data rows
        for i in 1:length(expt.τ_values)
            row = [i, expt.τ_values[i]]
            
            for region in expt.regions[]
                integrals = region.parameters[:integral].value[]
                errors = region.parameters[:error].value[]
                
                if i <= length(integrals) && i <= length(errors)
                    push!(row, integrals[i], errors[i])
                else
                    push!(row, "NaN", "NaN")
                end
            end
            
            println(f, join(row, "\t"))
        end
    end
end