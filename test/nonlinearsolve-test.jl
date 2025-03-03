using ForwardDiff
using NonlinearSolve
using LinearAlgebra
using Statistics
using GLMakie  # For visualization

include("nonlinearfit.jl")

# Generate some synthetic data with noise
function generate_data()
    a_true, b_true, c_true = 10.0, 0.5, 2.0
    x = collect(0:0.5:10)
    y_true = a_true .* exp.(-b_true .* x) .+ c_true
    y_noisy = y_true .+ 0.5 .* randn(length(x))
    return x, y_noisy, [a_true, b_true, c_true]
end

# Define our model and residual function
function model(x, p)
    a, b, c = p
    return a .* exp.(-b .* x) .+ c
end

# In-place residual function
function residual!(res, p, x, y)
    res .= model(x, p) .- y
end

# Main function to demonstrate the workflow
function fit_exponential_decay()
    # Generate synthetic data
    x, y, true_params = generate_data()
    
    # Initial parameter guess
    p0 = [8.0, 0.3, 1.0]  # [a, b, c]
    
    # Closure to capture x and y
    function resid!(res, p, _)
        residual!(res, p, x, y)
    end
    resid!(res, p) = resid!(res, p, nothing)
    
    # Pre-allocate a residual vector and determine its size
    resid_vec = similar(y)
    optimal_params, param_stderr = nonlinearfit(resid!, resid_vec, p0)
    
    # Print results
    println("True parameters: ", true_params)
    println("Optimal parameters: ", optimal_params)
    println("Parameter standard errors: ", param_stderr)
    
    # Calculate fitted curve and confidence intervals
    fitted_curve = model(x, optimal_params)
    
    # Create visualization with Makie
    fig = Figure()
    ax = Axis(fig[1, 1], 
              xlabel = "x", 
              ylabel = "y",
              title = "Exponential Decay Fit")
    
    # Plot the data points
    scatter!(ax, x, y, markersize = 8, color = :black, label = "Data")
    
    # Plot the fitted curve
    lines!(ax, x, fitted_curve, linewidth = 3, color = :blue, label = "Fitted curve")
    
    # Display the figure
    fig
    
    # return optimal_params, param_stderr, covariance, fig
end

# Run the example
fit_exponential_decay()

# Optional: save the figure
# save("exponential_fit.png", fig)