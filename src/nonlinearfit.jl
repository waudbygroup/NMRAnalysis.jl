"""
    leastsquresfit(resid!, resid_template, initial_params)

Perform nonlinear least squares fitting and calculate parameter uncertainties.

# Arguments
- `resid!(res, p)`: In-place residual function that modifies `res` based on parameters `p`
- `resid_template`: A pre-allocated residual vector to be overwritten
- `initial_params`: Initial guess for the parameters

# Returns
- `params`: Optimal parameter values
- `stderr`: Standard errors for each parameter or NaN values if covariance calculation fails
"""
function leastsquaresfit(resid!, resid_template, initial_params)
    # Set up problem for NonlinearSolve
    problem = NonlinearLeastSquaresProblem(
        NonlinearFunction(resid!, resid_prototype = resid_template), 
        initial_params
    )
    
    # Solve the problem
    sol = solve(problem)
    
    # Extract the best-fit parameters
    optimal_params = sol.u
    
    # Compute final residuals directly in the template
    resid!(resid_template, optimal_params)
    
    # Compute the Jacobian directly using ForwardDiff.jacobian
    jac = ForwardDiff.jacobian(resid!, resid_template, optimal_params)
    
    # Calculate statistical quantities
    sse = sum(resid_template.^2)  # Sum of squared errors
    dof = length(resid_template) - length(optimal_params)  # Degrees of freedom
    mse = sse / dof  # Mean squared error
    
    # Calculate covariance matrix with error trapping
    stderr = fill(0.0, length(optimal_params))
    try
        # Try using standard matrix inversion
        covariance = mse * inv(jac' * jac)
        stderr = sqrt.(diag(covariance))
    catch e
        # If standard inversion fails, try pseudoinverse
        try
            println("Warning: Standard matrix inversion failed. Attempting pseudoinverse.")
            covariance = mse * pinv(jac' * jac)
            stderr = sqrt.(diag(covariance))
        catch e2
            # If both methods fail, return NaN for stderr
            println("Warning: Could not calculate parameter uncertainties. Returning NaN values.")
            # stderr already filled with NaN values
        end
    end
    
    # Return only the optimal parameters and standard errors
    return optimal_params, stderr
end



"""
    curvefit(model, x, y, p0; weights=nothing, yerr=nothing)

Fit data to a nonlinear model using least squares optimization with optional weighting.

# Arguments
- `model`: A function of the form f(x, p) where x is the independent variable vector and p is the parameter vector
- `x`: Independent variable values (vector or matrix)
- `y`: Dependent variable measurements to fit
- `p0`: Initial parameter values for the model

# Keyword Arguments
- `weights`: Optional weights for weighted least squares. Cannot be used with `yerr`.
- `yerr`: Optional measurement errors in y. Cannot be used with `weights`. Internally converted to weights as 1/yerr^2.

# Returns
- `params`: Optimal parameter values
- `stderr`: Standard errors for each parameter or NaN values if covariance calculation fails
"""
function curvefit(model, x, y, p0; weights=nothing, yerr=nothing)
    # Check input dimensions
    if length(y) == 0
        error("Empty data array provided")
    end
    
    # Check that weights and yerr are not both provided
    if weights !== nothing && yerr !== nothing
        error("Cannot specify both `weights` and `yerr`. Choose only one weighting method.")
    end
    
    # Calculate weights from yerr if provided
    if yerr !== nothing
        # Check for zero or negative error values
        if any(yerr .<= 0)
            error("All values in `yerr` must be positive")
        end
        weights = 1 ./ (yerr.^2)
    end
    
    # Create residual function closure that captures the model, data, and weights
    function residual!(res, p)
        # Model accepts the entire x vector
        predictions = model(x, p)
        
        if weights !== nothing
            # For weighted fitting, scale the residuals by sqrt(weights)
            res .= sqrt.(weights) .* (y .- predictions)
        else
            res .= y .- predictions
        end
        
        return res
    end
    
    # Template for residuals with same size as y
    residual_template = similar(y)
    
    # Use the existing nonlinearfit function
    return leastsquresfit(residual!, residual_template, p0)
end