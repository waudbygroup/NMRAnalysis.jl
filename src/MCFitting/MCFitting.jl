module MCFitting

using MonteCarloMeasurements
using LsqFit
using Statistics
import MonteCarloMeasurements: ±

export fit_model, compare_models

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

        # Check if we need to recreate particles with a different count
        if length(y_data[1].particles) != n_particles
            # Create new particles with the requested particle count
            σ_vec = [pstd(y) for y in y_data]
            y_particles = [y_mean + σ_i * Particles(n_particles)
                           for (y_mean, σ_i) in zip(y_obs, σ_vec)]
        else
            y_particles = y_data
        end
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

    # fit to mean data first to get initial parameters
    # Use the mean of the particles for fitting
    y_mean = [pmean(y) for y in y_particles]
    fit_result = curve_fit(model_func, x_data, y_mean, initial_params)
    initial_params = coef(fit_result)

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

end