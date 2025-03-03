"""
    CustomModel <: ParametricModel

User-defined model with arbitrary function.
"""
struct CustomModel <: ParametricModel
    func::Function
    param_names::Vector{String}
    xlabel::String
    initial_values::Dict{String,Float64}
end

"""
    ExponentialModel <: ParametricModel

Exponential decay model: A * exp(-R * x)
"""
struct ExponentialModel <: ParametricModel
    func::Function
    param_names::Vector{String}
    xlabel::String
end

"""
    RecoveryModel <: ParametricModel

Recovery model: A * (1 - C * exp(-R * x))
"""
struct RecoveryModel <: ParametricModel
    func::Function
    param_names::Vector{String}
    xlabel::String
end

"""
    DiffusionModel <: ParametricModel

Diffusion model: A * exp(-(γ*δ*g)² * (Δ-δ/3) * D)
"""
struct DiffusionModel <: ParametricModel
    func::Function
    param_names::Vector{String}
    xlabel::String
    γ::Float64
    δ::Float64
    Δ::Float64
    σ::Float64
end

"""
    estimate_parameters(x::AbstractVector, y::AbstractVector, model::CustomModel)

Estimate initial parameters for custom model.
"""
function estimate_parameters(x::AbstractVector, y::AbstractVector, model::CustomModel)
    model.initial_values
end

"""
    ExponentialModel()

Create an exponential decay model: A * exp(-R * x)
"""
function ExponentialModel()
    ExponentialModel(
        (x, p) -> (@. p[1] * exp(-p[2] * x)),
        ["A", "R"],
        "Time / s"
    )
end

"""
    RecoveryModel()

Create a recovery model: A * (1 - C * exp(-R * x))
"""
function RecoveryModel()
    RecoveryModel(
        (x, p) -> (@. p[1] * (1 - p[2] * exp(-p[3] * x))),
        ["A", "C", "R"],
        "Time / s"
    )
end

"""
    DiffusionModel(γ, δ, Δ, σ=0.9)

Create a diffusion model: A * exp(-(γ*δ*σ*g)² * (Δ-δ/3) * D)
"""
function DiffusionModel(γ, δ, Δ, σ=0.9)
    DiffusionModel(
        (g, p) -> (@. p[1] * exp(-(γ * δ * σ * g)^2 * (Δ - δ/3) * p[2])),
        ["A", "D"],
        "Gradient strength / T m⁻¹",
        γ, δ, Δ, σ
    )
end

"""
    CustomModel(modelfunction::String, params::Vector{Pair{String,Float64}}, xlabel="x")

Create a custom model with user-defined function and parameters.
"""
function CustomModel(modelfunction::String, params::Vector{Pair{String,Float64}}, xlabel="x")
    param_names = first.(params)
    CustomModel(
        modeltofunction(modelfunction, param_names),
        param_names,
        xlabel,
        Dict(params)
    )
end

"""
    modeltofunction(expr, param_names)

Convert a string expression to a function with parameters.
"""
function modeltofunction(expr, param_names)
    function replace_walker(ex)
        if ex isa Symbol && string(ex) in param_names
            # Found a parameter name as a symbol - replace with indexed access
            idx = findfirst(==(string(ex)), param_names)
            return :(p[$idx])
        elseif ex isa Expr
            # Recursively walk through expression
            return Expr(ex.head, map(replace_walker, ex.args)...)
        else
            # Leave other elements (numbers, etc) unchanged
            return ex
        end
    end
    parsed_expr = Meta.parse(expr)
    expr_with_params = replace_walker(parsed_expr)
    func_expr = :($(Expr(:tuple, :x, :p)) -> @. $expr_with_params)
    eval(func_expr)
end

"""
    postfit!(region::Region, expt::Experiment1D, ::NoFitting)

Post-fit processing for no-fitting model.
"""
function postfit!(region::Region, expt::Experiment1D, ::NoFitting)
    region.postfitted[] = true
end

"""
    postfit!(region::Region, expt::Experiment1D, model::ParametricModel)

Generic post-fit processing for parametric models.
"""
function postfit!(region::Region, expt::Experiment1D, model::ParametricModel)
    @debug "Post-fitting model"
    
    # Get x and y data
    # (Each experiment type should provide the appropriate data)
    x, y = get_fit_data(region, expt)
    
    # Initial parameter estimates
    p0 = collect(values(estimate_parameters(x, y, model)))
    
    # Fit the model
    fit = curve_fit(model.func, x, y, p0)
    pfit = coef(fit)
    perr = stderror(fit)
    
    # Update post-parameters with fitted values
    for (i, name) in enumerate(model.param_names)
        param = region.postparameters[Symbol(name)]
        param.value[] .= pfit[i]
        param.uncertainty[] .= perr[i]
    end
    
    @debug "Fitted parameters: $(region.postparameters)"
    region.postfitted[] = true
end

"""
    get_model_data(region, expt::Experiment1D, ::NoFitting)

Get data for visualization with no fitting model.
"""
function get_model_data(region, expt::Experiment1D, ::NoFitting)
    isnothing(region) && return (Point2f[], [(0.0, 0.0, 0.0)], Point2f[])
    
    x, y = get_fit_data(region, expt)
    
    # Create points for plot
    obs_points = Point2f.(x, y)
    obs_errors = [(x[i], y[i], 0.0) for i in 1:length(y)]  # No error estimates
    fit_points = Point2f[]  # No fit line
    
    return (obs_points, obs_errors, fit_points)
end

"""
    get_model_data(region, expt::Experiment1D, model::ParametricModel)

Get data for visualization with parametric model.
"""
function get_model_data(region, expt::Experiment1D, model::ParametricModel)
    isnothing(region) && return (Point2f[], [(0.0, 0.0, 0.0)], Point2f[])
    
    x, y = get_fit_data(region, expt)
    
    # Create points for plot
    obs_points = Point2f.(x, y)
    
    # Use measurement errors if available
    if haskey(region.parameters, :error)
        err = region.parameters[:error].value[]
        obs_errors = [(x[i], y[i], err[i]) for i in 1:length(y)]
    else
        obs_errors = [(x[i], y[i], 0.0) for i in 1:length(y)]
    end
    
    # Calculate fit line if region has been fitted
    if region.postfitted[]
        xpred = range(min(0.0, minimum(x)), 1.1*maximum(x), 100)
        p = [region.postparameters[Symbol(name)].value[][1] for name in model.param_names]
        ypred = model.func(xpred, p)
        fit_points = Point2f.(xpred, ypred)
    else
        fit_points = Point2f[]
    end
    
    return (obs_points, obs_errors, fit_points)
end

"""
    model_parameter_text(region::Region, ::NoFitting)

Get text describing parameters for no-fitting model.
"""
function model_parameter_text(region::Region, ::NoFitting)
    if haskey(region.parameters, :integral)
        integral = region.parameters[:integral].value[]
        if integral isa AbstractArray && length(integral) > 1
            return ["Integral: $(integral[1]) (slice 1)"]
        else
            return ["Integral: $integral"]
        end
    else
        return ["No integral calculated"]
    end
end

"""
    model_parameter_text(region::Region, model::ParametricModel)

Get text describing parameters for parametric model.
"""
function model_parameter_text(region::Region, model::ParametricModel)
    map(model.param_names) do name
        param = region.postparameters[Symbol(name)]
        "$name: $(param.value[][1]) ± $(param.uncertainty[][1])"
    end
end

"""
    model_info_text(::NoFitting, x::AbstractVector)

Get text describing model info for no-fitting model.
"""
function model_info_text(::NoFitting, x::AbstractVector)
    ["Number of spectra: $(length(x))"]
end

"""
    model_info_text(model::ParametricModel, x::AbstractVector)

Get text describing model info for parametric model.
"""
function model_info_text(model::ParametricModel, x::AbstractVector)
    [
        "Number of points: $(length(x))",
        "$(model.xlabel) range: $(minimum(x)) - $(maximum(x))"
    ]
end

"""
    get_fit_data(region::Region, expt::Experiment1D)

Get x and y data for fitting. Should be implemented by each experiment type.
"""
function get_fit_data(region::Region, expt::Experiment1D)
    # Default implementation returns an empty dataset
    # Each experiment type should override this
    return Float64[], Float64[]
end

"""x::AbstractVector, y::AbstractVector, ::RecoveryModel)

Estimate initial parameters for recovery model.
"""
function estimate_parameters(x::AbstractVector, y::AbstractVector, ::RecoveryModel)
    A = maximum(y)
    C = 1.0
    R = 3.0 / maximum(x)
    Dict("A" => A, "C" => C, "R" => R)
end

"""
    estimate_parameters(x::AbstractVector, y::AbstractVector, ::ExponentialModel)

Estimate initial parameters for exponential decay model.
"""
function estimate_parameters(x::AbstractVector, y::AbstractVector, ::ExponentialModel)
    A = maximum(abs.(y))
    R = 3.0 / maximum(x)
    Dict("A" => A, "R" => R)
end

"""
    estimate_parameters(x::AbstractVector, y::AbstractVector, model::DiffusionModel)

Estimate initial parameters for diffusion model.
"""
function estimate_parameters(x::AbstractVector, y::AbstractVector, ::DiffusionModel)
    A = maximum(abs.(y))
    D = 1.0e-10  # Initial guess for diffusion coefficient (m²/s)
    Dict("A" => A, "D" => D)
end

