# Functions for model fitting with MC uncertainty and model comparison
using MonteCarloMeasurements
using LsqFit
using Statistics
import MonteCarloMeasurements: ±

±(μ, σ, n=2000) = μ + σ * repeat(Particles(n), length(μ))  # Create a particle distribution

"""
    fit_model(x_data, y_data, model_func, initial_params; 
                              σ=nothing, n_particles=2000, param_names=nothing)

Fit a model to data with Monte Carlo uncertainty propagation.

# Arguments
- `x_data`: Independent variable values
- `y_data`: Dependent variable values (can be raw observations or Particles with uncertainty)
- `model_func`: Function that takes x and parameters as arguments
- `initial_params`: Initial parameter values for fitting
- `σ`: Optional standard deviation(s) of observations (scalar or vector) - not needed if y_data has uncertainty
- `n_particles`: Number of particles for uncertainty propagation (only used with raw observations)
- `param_names`: Optional names for parameters (for better reporting)

# Returns
- Named tuple with fitted parameters, predicted values, log-likelihood, AIC, etc.
"""
function fit_model(x_data, y_data, model_func, initial_params;
                   σ=nothing, n_particles=2000, param_names=nothing)
    # Determine if y_data already has uncertainty
    has_uncertainty = eltype(y_data) <: Particles

    # Extract mean values for reporting and log-likelihood calculation
    if has_uncertainty
        y_obs = [pmean(y) for y in y_data]
        y_particles = y_data
    else
        y_obs = y_data
        # Create particle observations if needed
        if isnothing(σ)
            error("Either provide y_data with uncertainty or specify σ")
        elseif isa(σ, Number)
            y_particles = y_obs .+ σ * [Particles(n_particles) for _ in 1:length(y_obs)]
            σ_vec = fill(σ, length(y_obs))
        else
            y_particles = [y + σ_i * Particles(n_particles) for (y, σ_i) in zip(y_obs, σ)]
            σ_vec = σ
        end
    end

    # Extract σ values for log-likelihood if y_data has uncertainty
    if has_uncertainty && isnothing(σ)
        σ_vec = [pstd(y) for y in y_data]
    end

    # Create fit function (closure)
    function fit_from_observations(y_data)
        fit_result = curve_fit(model_func, x_data, y_data, initial_params)
        return coef(fit_result)
    end

    # Fit the model to all particles
    fitted_params_raw = @bymap fit_from_observations(y_particles)

    # Calculate predicted values
    predicted_values = model_func(x_data, fitted_params_raw)

    # Calculate log-likelihood and AIC
    n = length(y_obs)
    k = length(initial_params)

    # Calculate log-likelihood for each particle
    log_likelihood = -n / 2 * log(2π) - sum(log.(σ_vec)) -
                     0.5 *
                     sum([(y_obs[i] - predicted_values[i])^2 / σ_vec[i]^2 for i in 1:n])

    # Calculate AIC
    aic = 2 * k - 2 * log_likelihood

    # Parameter names for reporting
    if isnothing(param_names)
        param_names = ["p$i" for i in 1:length(initial_params)]
    end

    # Prepare parameter summary
    param_summary = [(name=param_names[i],
                      value=fitted_params_raw[i],
                      mean=pmean(fitted_params_raw[i]),
                      std=pstd(fitted_params_raw[i])) for i in 1:length(fitted_params_raw)]

    return (parameters=fitted_params_raw,
            param_summary=param_summary,
            predicted=predicted_values,
            residuals=y_obs .- [pmean(y) for y in predicted_values],
            log_likelihood=log_likelihood,
            aic=aic,
            n_params=k,
            data=(x=x_data, y=y_obs, σ=σ_vec))
end

"""
    compare_models(fit_result1, fit_result2; confidence_threshold=0.95)

Compare two model fits using AIC and uncertainty propagation.

# Arguments
- `fit_result1`: Result from first model fit
- `fit_result2`: Result from second model fit
- `confidence_threshold`: Confidence threshold for model rejection

# Returns
- Named tuple with comparison results
"""
function compare_models(fit_result1, fit_result2; confidence_threshold=0.95)
    # Calculate ΔAIC (positive means model2 is better)
    delta_aic = fit_result1.aic - fit_result2.aic

    # Calculate probability that model2 is better
    prob_model2_better = mean(delta_aic.particles .> 0)

    # Determine if we can reject model1 in favor of model2
    reject_model1 = prob_model2_better > confidence_threshold

    # Determine if we can reject model2 in favor of model1
    reject_model2 = mean(delta_aic.particles .< 0) > confidence_threshold

    # Determine if models are different in complexity
    n_params1, n_params2 = fit_result1.n_params, fit_result2.n_params

    # Check if simpler model is rejected in favor of more complex
    simpler_rejected = false
    if n_params1 < n_params2
        simpler_rejected = reject_model1
    elseif n_params2 < n_params1
        simpler_rejected = reject_model2
    end

    return (delta_aic=delta_aic,
            prob_model2_better=prob_model2_better,
            reject_model1=reject_model1,
            reject_model2=reject_model2,
            simpler_rejected=simpler_rejected,
            avg_delta_aic=pmean(delta_aic),
            std_delta_aic=pstd(delta_aic))
end

"""
    summarize_fit(fit_result)

Generate a text summary of the model fit results.

# Arguments
- `fit_result`: Result from model fit

# Returns
- String with summary information
"""
function summarize_fit(fit_result)
    io = IOBuffer()
    println(io, "Model Fit Summary")
    println(io, "=================")
    println(io, "Number of parameters: $(fit_result.n_params)")
    println(io, "AIC: $(pmean(fit_result.aic)) ± $(pstd(fit_result.aic))")
    println(io,
            "Log-likelihood: $(pmean(fit_result.log_likelihood)) ± $(pstd(fit_result.log_likelihood))")

    println(io, "\nParameters:")
    for p in fit_result.param_summary
        println(io,
                "  $(rpad(p.name, 15)) = $(round(p.mean, digits=4)) ± $(round(p.std, digits=4))")
    end

    println(io, "\nResidual statistics:")
    println(io,
            "  Mean absolute residual: $(round(mean(abs.(fit_result.residuals)), digits=6))")
    println(io,
            "  Root mean square residual: $(round(sqrt(mean(fit_result.residuals.^2)), digits=6))")

    return String(take!(io))
end

"""
    summarize_comparison(comparison_result, model1_name, model2_name)

Generate a text summary of the model comparison results.

# Arguments
- `comparison_result`: Result from model comparison
- `model1_name`: Name of the first model
- `model2_name`: Name of the second model

# Returns
- String with summary information
"""
function summarize_comparison(comparison_result, model1_name, model2_name)
    io = IOBuffer()
    println(io, "Model Comparison: $model1_name vs $model2_name")
    println(io, "==========================================")

    delta = comparison_result.avg_delta_aic
    if delta > 0
        better_model = model2_name
        worse_model = model1_name
    else
        better_model = model1_name
        worse_model = model2_name
        delta = -delta
    end

    println(io,
            "ΔAIC: $(round(abs(delta), digits=2)) ± $(round(comparison_result.std_delta_aic, digits=2))")
    println(io,
            "Probability that $model2_name is better: $(round(comparison_result.prob_model2_better, digits=4))")

    if comparison_result.reject_model1
        println(io,
                "\nConclusion: $model1_name can be rejected in favor of $model2_name (p = $(round(comparison_result.prob_model2_better, digits=4)))")
    elseif comparison_result.reject_model2
        println(io,
                "\nConclusion: $model2_name can be rejected in favor of $model1_name (p = $(round(1 - comparison_result.prob_model2_better, digits=4)))")
    else
        println(io, "\nConclusion: Cannot confidently choose between models")
    end

    if comparison_result.simpler_rejected
        println(io, "The simpler model is rejected in favor of the more complex model.")
    end

    # AIC interpretation guide
    println(io, "\nAIC Interpretation Guide:")
    println(io, "  ΔAIC < 2: Models are essentially equivalent")
    println(io, "  ΔAIC 2-4: There is some evidence for the better model")
    println(io, "  ΔAIC 4-7: There is considerable evidence for the better model")
    println(io, "  ΔAIC > 10: There is strong evidence for the better model")

    return String(take!(io))
end

"""
    plot_fit_comparison(fit_result1, fit_result2, model1_name, model2_name; 
                       true_y=nothing, size=(900, 600))

Plot a comparison of two model fits.

# Arguments
- `fit_result1`: Result from first model fit
- `fit_result2`: Result from second model fit
- `model1_name`: Name of the first model
- `model2_name`: Name of the second model
- `true_y`: Optional true y values for comparison
- `size`: Size of the figure

# Returns
- Figure object
"""
function plot_fit_comparison(fit_result1, fit_result2, model1_name, model2_name;
                             true_y=nothing, size=(900, 600))

    # Extract data
    x_data = fit_result1.data.x
    y_obs = fit_result1.data.y
    σ = fit_result1.data.σ

    # Create figure
    fig = Figure(; size=size)

    # Plot fits
    ax1 = Axis(fig[1, 1:2]; title="Model Comparison", xlabel="x", ylabel="y")

    # Plot observed data
    scatter!(ax1, x_data, y_obs; color=:black, markersize=6, label="Observed")

    # Plot error bars if σ is a vector
    if isa(σ, AbstractVector)
        errorbars!(ax1, x_data, y_obs, σ; color=:black, whiskerwidth=10)
    end

    # Plot model 1 fit
    y_pred1 = [pmean(y) for y in fit_result1.predicted]
    lower1 = [pmean(y) - 2 * pstd(y) for y in fit_result1.predicted]
    upper1 = [pmean(y) + 2 * pstd(y) for y in fit_result1.predicted]
    band!(ax1, x_data, lower1, upper1; color=(:red, 0.2))
    lines!(ax1, x_data, y_pred1; color=:red, label=model1_name)

    # Plot model 2 fit
    y_pred2 = [pmean(y) for y in fit_result2.predicted]
    lower2 = [pmean(y) - 2 * pstd(y) for y in fit_result2.predicted]
    upper2 = [pmean(y) + 2 * pstd(y) for y in fit_result2.predicted]
    band!(ax1, x_data, lower2, upper2; color=(:blue, 0.2))
    lines!(ax1, x_data, y_pred2; color=:blue, label=model2_name)

    # Plot true y if provided
    if !isnothing(true_y)
        lines!(ax1, x_data, true_y; color=:black, linestyle=:dash, label="True")
    end

    axislegend(ax1; position=:rt)

    # Plot residuals
    ax2 = Axis(fig[2, 1]; title="Residuals (Model 1)", xlabel="x", ylabel="Residual")
    scatter!(ax2, x_data, fit_result1.residuals; color=:red)
    hlines!(ax2, [0]; color=:black, linestyle=:dash)

    ax3 = Axis(fig[2, 2]; title="Residuals (Model 2)", xlabel="x", ylabel="Residual")
    scatter!(ax3, x_data, fit_result2.residuals; color=:blue)
    hlines!(ax3, [0]; color=:black, linestyle=:dash)

    # Link axes
    linkyaxes!(ax2, ax3)

    return fig
end

"""
    plot_parameter_distributions(fit_result; size=(900, 600), 
                                true_params=nothing)

Plot parameter distributions from a model fit.

# Arguments
- `fit_result`: Result from model fit
- `size`: Size of the figure
- `true_params`: Optional true parameter values

# Returns
- Figure object
"""
function plot_parameter_distributions(fit_result; size=(900, 600),
                                      true_params=nothing)
    # Create figure
    fig = Figure(; size=size)

    # Get number of parameters
    n_params = length(fit_result.parameters)

    # Calculate grid layout
    n_cols = min(3, n_params)
    n_rows = ceil(Int, n_params / n_cols)

    # Create parameter plots
    for i in 1:n_params
        row = (i - 1) ÷ n_cols + 1
        col = (i - 1) % n_cols + 1

        ax = Axis(fig[row, col];
                  title=fit_result.param_summary[i].name,
                  xlabel="Value", ylabel="Density")

        hist!(ax, fit_result.parameters[i].particles;
              bins=30, color=(:blue, 0.5),
              normalization=:pdf)

        # Add mean line
        vlines!(ax, [fit_result.param_summary[i].mean];
                color=:red, linewidth=2,
                label="Mean")

        # Add true parameter if provided
        if !isnothing(true_params) && length(true_params) >= i
            vlines!(ax, [true_params[i]];
                    color=:black, linestyle=:dash,
                    label="True")
        end

        if i == 1
            axislegend(ax; position=:rt)
        end
    end

    return fig
end

"""
    plot_aic_comparison(fit_result1, fit_result2, model1_name, model2_name; 
                       size=(900, 600))

Plot AIC distributions for two model fits.

# Arguments
- `fit_result1`: Result from first model fit
- `fit_result2`: Result from second model fit
- `model1_name`: Name of the first model
- `model2_name`: Name of the second model
- `size`: Size of the figure

# Returns
- Figure object
"""
function plot_aic_comparison(fit_result1, fit_result2, model1_name, model2_name;
                             size=(900, 600))
    # Calculate ΔAIC
    delta_aic = fit_result1.aic - fit_result2.aic
    prob = mean(delta_aic.particles .> 0)

    # Create figure
    fig = Figure(; size=size)

    # Plot AIC distributions
    ax1 = Axis(fig[1, 1]; title="AIC Distributions",
               xlabel="AIC", ylabel="Density")

    hist!(ax1, fit_result1.aic.particles; bins=30,
          color=(:red, 0.5), normalization=:pdf,
          label=model1_name)

    hist!(ax1, fit_result2.aic.particles; bins=30,
          color=(:blue, 0.5), normalization=:pdf,
          label=model2_name)

    axislegend(ax1; position=:rt)

    # Plot ΔAIC distribution
    ax2 = Axis(fig[1, 2];
               title="ΔAIC Distribution ($model1_name - $model2_name)",
               xlabel="ΔAIC", ylabel="Density")

    hist!(ax2, delta_aic.particles; bins=30,
          color=(:purple, 0.5), normalization=:pdf)

    vlines!(ax2, [0]; color=:black, linestyle=:dash)
    vlines!(ax2, [pmean(delta_aic)]; color=:red,
            linewidth=2, label="Mean")

    # Add text annotation for probability
    text!(ax2, "P(ΔAIC > 0) = $(round(prob, digits=4))";
          position=(0.7, 0.9), align=(:right, :center),
          space=:relative, fontsize=16)

    axislegend(ax2; position=:rt)

    return fig
end