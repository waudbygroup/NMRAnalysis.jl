using MonteCarloMeasurements
using LsqFit
using GLMakie
import MonteCarloMeasurements: ±

±(μ, σ) = μ + σ * repeat(Particles(2000), length(μ))  # Create a particle distribution

function create_fit_function(x_data, model_func, initial_params)
    # Return a closure that only needs y_data as input
    function fit_from_observations(y_data)
        fit_result = curve_fit(model_func, x_data, y_data, initial_params)
        return coef(fit_result)  # Return just the fitted parameters
        # Alternative: return (parameters = coef(fit_result), errors = stderror(fit_result))
    end

    return fit_from_observations
end

# Example usage:
function exponential_decay(x, p)
    return @. p[1] * exp(-p[2] * x)
end

x = 0:0.5:5
ptrue = [2.0, 0.5]
ytrue = exponential_decay(x, ptrue)

# Create the fit function
p0 = [1.0, 0.1]
σ = 0.2
fit_func = create_fit_function(x, exponential_decay, p0)

# Now this can be used with MonteCarloMeasurements
yobs = ytrue .+ randn(length(x)) * σ  # Add noise
yparticles = ±(yobs, σ)  # Create a particle distribution

# Map the fit function across particle realizations
fitresults = @bymap fit_func(yparticles)
Afit = fitresults[1]
kfit = fitresults[2]
As = Afit.particles
ks = kfit.particles
@info "Fitted parameters" Afit kfit

yfit = exponential_decay(x, fitresults)  # Fitted curve

# # Calculate the log-likelihood and AIC
# n = length(x)
# k = length(p0)
# RSS = sum((yobs .- yfit).^2)
# LL = -n/2 * log(2π*σ^2) - RSS/(2σ^2)
# AIC = 2*k - 2*LL
# @info "LL and AIC" LL AIC

# plot results
fig = Figure(; size=(1000, 500))
ax1 = Axis(fig[1, 1]; title="Fitted Curve", xlabel="x", ylabel="y")
lines!(ax1, x, ytrue; color=:dodgerblue, label="True")
band!(ax1, x, yfit; color=:red, alpha=0.2, label="Fit")
rangebars!(ax1, x, yparticles; color=:red, whiskerwidth=10)
scatter!(ax1, x, yobs; color=:red, label="Observed")
axislegend(ax1; position=:rt)

# plot 2d histogram of the fitted parameters
ax2 = Axis(fig[1, 2]; title="Fitted Parameters", xlabel="A", ylabel="k")
hexbin!(ax2, As, ks; colormap=:viridis)
Amean = pmean(Afit)
kmean = pmean(kfit)
Astd = pstd(Afit)
kstd = pstd(kfit)
errorbars!(ax2, [Amean], [kmean], [Astd]; direction=:x, color=:orangered, whiskerwidth=10)
errorbars!(ax2, [Amean], [kmean], [kstd]; direction=:y, color=:orangered, whiskerwidth=10)
scatter!(ax2, Point2f(Amean, kmean); color=:orangered, label="Fit")
scatter!(ax2, Point2f(ptrue...); color=:magenta, marker=:star5, markersize=15, label="True")

# # plot LL distribution
# ax3 = Axis(fig[1,3], title = "Log-Likelihood Distribution", xlabel= "Log-Likelihood", ylabel = "Density")
# hist!(ax3, LL.particles)

# return the figure
fig
