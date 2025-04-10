using MonteCarloMeasurements
using LsqFit
using GLMakie
using Statistics

# Define model functions
function exponential_decay(x, p)
    return @. p[1] * exp(-exp(p[2]) * x)
end

function constant_model(x, p)
    return fill(p[1], length(x))
end

# Example 1: Using raw data with uncertainty parameter
function example_with_raw_data()
    println("\n=== Example 1: Using raw data with uncertainty parameter ===\n")

    # Generate synthetic data
    x = 0:0.5:5
    true_params = [1.20, -0.5]
    ytrue = exponential_decay(x, true_params)
    noise_level = 0.1
    yobs = ytrue .+ randn(length(x)) * noise_level

    # Fit models
    exp_fit = fit_model(x, yobs,
                        exponential_decay,
                        [1.0, 0.0];  # Initial parameters 
                        σ=noise_level,
                        param_names=["Amplitude", "log(Rate)"])

    const_fit = fit_model(x, yobs,
                          constant_model,
                          [0.5];  # Initial parameters
                          σ=noise_level,
                          param_names=["Constant"])

    # Compare models
    comparison = compare_models(const_fit, exp_fit)

    # Print summaries
    println("Exponential Model Fit:")
    println(summarize_fit(exp_fit))
    println("\nConstant Model Fit:")
    println(summarize_fit(const_fit))
    println("\nModel Comparison:")
    println(summarize_comparison(comparison, "Constant", "Exponential"))

    # Create plots
    fig_fits = plot_fit_comparison(const_fit, exp_fit,
                                   "Constant", "Exponential";
                                   true_y=ytrue)

    display(fig_fits)
    readline()

    fig_params = plot_parameter_distributions(exp_fit;
                                              true_params=true_params)

    display(fig_params)
    readline()

    fig_aic = plot_aic_comparison(const_fit, exp_fit,
                                  "Constant", "Exponential")

    return display(fig_aic)
end

# Example 2: Using data with uncertainty already included
function example_with_particle_data()
    println("\n=== Example 2: Using data with uncertainty already included ===\n")

    # Generate synthetic data with uncertainty
    x = 0:0.5:5
    true_params = [1.20, -0.5]
    ytrue = exponential_decay(x, true_params)
    noise_level = 0.1

    # Create particle distributions directly
    yobs = ytrue .+ noise_level * randn(length(ytrue))
    y_particles = ±(yobs, noise_level)

    # Fit models
    exp_fit = fit_model(x, y_particles,
                        exponential_decay,
                        [1.0, 0.0];  # Initial parameters
                        # σ=σ_values,  # Still need σ for likelihood calculation
                        param_names=["Amplitude", "log(Rate)"])

    const_fit = fit_model(x, y_particles,
                          constant_model,
                          [0.5];  # Initial parameters
                          # σ=σ_values,  # Still need σ for likelihood calculation
                          param_names=["Constant"])

    # Compare models
    comparison = compare_models(const_fit, exp_fit)

    # Print summaries
    println("Exponential Model Fit:")
    println(summarize_fit(exp_fit))
    println("\nConstant Model Fit:")
    println(summarize_fit(const_fit))
    println("\nModel Comparison:")
    return println(summarize_comparison(comparison, "Constant", "Exponential"))
end

# Run the examples
example_with_raw_data()
example_with_particle_data()