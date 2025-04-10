using MonteCarloMeasurements
using LsqFit
using GLMakie
using Statistics
import MonteCarloMeasurements: ±

±(μ, σ) = μ + σ * repeat(Particles(2000), length(μ))  # Create a particle distribution

function create_fit_function(x_data, model_func, initial_params)
    # Return a closure that only needs y_data as input
    function fit_from_observations(y_data)
        fit_result = curve_fit(model_func, x_data, y_data, initial_params)
        return coef(fit_result)  # Return just the fitted parameters
    end

    return fit_from_observations
end

# Define models
function exponential_decay(x, p)
    return @. p[1] * exp(-exp(p[2]) * x)
end

function constant_model(x, p)
    return fill(p[1], length(x))
end

# Generate data
x = 0:0.5:5
ptrue = [1.20, -0.5]
ytrue = exponential_decay(x, ptrue)

# Set up fitting
σ = 0.1
nparticles = 100
yobs = ytrue .+ randn(length(x)) * σ  # Add noise
yparticles = yobs + σ * [Particles(nparticles) for i in 1:length(yobs)]
# Exponential model fitting
p0_exp = [1.0, 0]
fit_func_exp = create_fit_function(x, exponential_decay, p0_exp)
fitresults_exp = @bymap fit_func_exp(yparticles)
Afit = fitresults_exp[1]
kfit = fitresults_exp[2]
As = Afit.particles
ks = kfit.particles
yfit_exp = exponential_decay(x, fitresults_exp)

# Constant model fitting
p0_const = [mean(yobs)]
fit_func_const = create_fit_function(x, constant_model, p0_const)
fitresults_const = @bymap fit_func_const(yparticles)
cfit = fitresults_const[1]
cs = cfit.particles
yfit_const = constant_model(x, fitresults_const)

# Calculate log-likelihood for each ensemble member
n = length(yobs)
LL_exp = -n / 2 * log(2π * σ^2) - sum((yobs .- yfit_exp) .^ 2) / (2σ^2)
LL_const = -n / 2 * log(2π * σ^2) - sum((yobs .- yfit_const) .^ 2) / (2σ^2)

# Calculate AIC for both models
k_exp = length(p0_exp)  # Number of parameters in exponential model
k_const = length(p0_const)  # Number of parameters in constant model

AIC_exp = 2 * k_exp - 2 * LL_exp
AIC_const = 2 * k_const - 2 * LL_const

# Calculate differences in AIC (positive means exponential is better)
ΔAIC = AIC_const - AIC_exp
p_exp_better = mean(ΔAIC.particles .> 0)
reject_const = p_exp_better > 0.95

# Print results
@info "Model comparison" ΔAIC p_exp_better
@info "Reject constant model" reject_const

# Create plots
fig = Figure(; size=(1200, 800))

# Plot fitted curves
ax1 = Axis(fig[1, 1:2]; title="Fitted Curves", xlabel="x", ylabel="y")
lines!(ax1, x, ytrue; color=:dodgerblue, label="True")
band!(ax1, x, yfit_exp; color=:red, alpha=0.2, label="Exponential fit")
band!(ax1, x, yfit_const; color=:green, alpha=0.2, label="Constant fit")
scatter!(ax1, x, yobs; color=:black, label="Observed")
rangebars!(ax1, x, yparticles; color=:black, whiskerwidth=10)
axislegend(ax1; position=:rt)

# Plot parameter distribution for exponential model
ax2 = Axis(fig[1, 3]; title="Exponential Model Parameters", xlabel="A", ylabel="k")
hexbin!(ax2, As, ks; colormap=:viridis)
Amean = pmean(Afit)
kmean = pmean(kfit)
Astd = pstd(Afit)
kstd = pstd(kfit)
errorbars!(ax2, [Amean], [kmean], [Astd]; direction=:x, color=:orangered, whiskerwidth=10)
errorbars!(ax2, [Amean], [kmean], [kstd]; direction=:y, color=:orangered, whiskerwidth=10)
scatter!(ax2, Point2f(Amean, kmean); color=:orangered, label="Fit")
scatter!(ax2, Point2f(ptrue...); color=:magenta, marker=:star5, markersize=15, label="True")

# Plot log-likelihood distributions
ax3 = Axis(fig[2, 1]; title="Log-Likelihood Distributions", xlabel="Log-Likelihood",
           ylabel="Density")
hist!(ax3, all_LL_exp; bins=30, color=(:red, 0.5), label="Exponential")
hist!(ax3, all_LL_const; bins=30, color=(:green, 0.5), label="Constant")
axislegend(ax3; position=:rt)

# Plot AIC distributions
ax4 = Axis(fig[2, 2]; title="AIC Distributions", xlabel="AIC", ylabel="Density")
hist!(ax4, all_AIC_exp; bins=30, color=(:red, 0.5), label="Exponential")
hist!(ax4, all_AIC_const; bins=30, color=(:green, 0.5), label="Constant")
axislegend(ax4; position=:rt)

# Plot ΔAIC distribution
ax5 = Axis(fig[2, 3]; title="ΔAIC Distribution (Constant - Exponential)", xlabel="ΔAIC",
           ylabel="Density")
hist!(ax5, ΔAIC.particles; bins=30)
vlines!(ax5, [0]; color=:black, linestyle=:dash)

fig