struct CustomModel <: ParametricModel
    func::Function
    param_names::Vector{String}
    xlabel::String
    initial_values::Dict{String,Float64}
end

struct ExponentialModel <: ParametricModel
    func::Function
    param_names::Vector{String}
    xlabel::String
end

struct RecoveryModel <: ParametricModel
    func::Function
    param_names::Vector{String}
    xlabel::String
end

struct DiffusionModel <: ParametricModel
    func::Function
    param_names::Vector{String}
    xlabel::String
    q::Float64
end

function estimate_parameters(x, y, ::RecoveryModel)
    A = maximum(y)
    C = 1.0
    R = 3.0 / maximum(x)
    Dict("A" => A, "C" => C, "R" => R)
end

function estimate_parameters(x, y, ::ExponentialModel)
    A = maximum(abs.(y))
    R = 3.0 / maximum(x)  # Different heuristic for T1
    Dict("A" => A, "R" => R)
end

function estimate_parameters(_, _, model::CustomModel)
    model.initial_values
end

function ExponentialModel()
    ExponentialModel(
        (x, p) -> (@. p[1] * exp(-p[2] * x)),
        ["A", "R"],
        "Time / s"
    )
end

function RecoveryModel()
    RecoveryModel(
        (x, p) -> (@. p[1] * (1 - p[2] * exp(-p[3] * x))),
        ["A", "C", "R"],
        "Time / s"
    )
end

function DiffusionModel(δ, Δ, σ, γ)
    q = (γ*δ*σ)^2 * (Δ - δ/3)
    DiffusionModel(
        (g, p) -> (@. p[1] * exp(-q * 1e10 * p[2] * g^2)),
        ["A", "D_e10", "R"],
        "Gradient strength / T m⁻¹",
        q
    )
end

function CustomModel(modelfunction::String, params::Vector{Pair{String,Float64}}, xlabel="x")
    param_names = first.(params)
    CustomModel(
        modeltofunction(modelfunction, param_names),
        param_names,
        xlabel,
        Dict(params)
    )
end

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

function modelinfo(::NoFitting, x)
    ["Number of spectra: $(length(x))"]
end

function modelinfo(model::ParametricModel, x)
    [
        "Number of spectra: $(length(x))",
        "$(model.xlabel) range: $(minimum(x)) - $(maximum(x))"
    ]
end

function modelparameterinfo(_::Region, ::NoFitting)
    []
end

function modelparameterinfo(region::Region, model::ParametricModel)
    map(model.param_names) do name
        param = region.fitparameters[Symbol(name)]
        "$name: $(param.value[][1] ± param.uncertainty[][1])"
    end
end
