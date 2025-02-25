model_R1rho_onres(νSL, p) = @. p[2] + p[3] * exp(2p[4]) / (exp(2p[4]) + (2π * νSL)^2)
model_I_onres(TSL, νSL, p) = p[1] * exp.(-TSL .* model_R1rho_onres(νSL, p))
model_I_onres(x, p) = model_I_onres(x[:,1], x[:,2], p)

function fit_onres(state, p0)
    dataset = state[:dataset]

    x = hcat(dataset.TSLs, dataset.νSLs)
    intensities = state[:intensities][]
    
    @debug "running fit..." p0
    fit = LsqFit.curve_fit(model_I_onres, x, intensities, p0)
    @debug fit.param fit.converged
    return fit
end

# Null model for R1rho (no exchange)
model_R1rho_onres_null(νSL, p) = @. p[2]  # Only R2,0 term, no Rex contribution
model_I_onres_null(TSL, νSL, p) = p[1] * exp.(-TSL .* model_R1rho_onres_null(νSL, p))
model_I_onres_null(x, p) = model_I_onres_null(x[:,1], x[:,2], p)

# p0 = [I0, R2,0]
function fit_onres_null(state, p0)
    dataset = state[:dataset]

    x = hcat(dataset.TSLs, dataset.νSLs)
    intensities = state[:intensities][]
    
    @debug "running null model fit..." p0
    # Only need parameters p[1] (I0) and p[2] (R2,0)
    fit = LsqFit.curve_fit(model_I_onres_null, x, intensities, p0)
    @debug fit.param fit.converged
    return fit
end

function fitexp(t, I, p)
    I0 = p[1]
    R0 = p[3] / 2
    p0 = [I0, R0]

    model(t, p) = @. p[1] * exp(-t * p[2])

    R, Re = try
        fit = LsqFit.curve_fit(model, t, I, p0)
        R = coef(fit)[2]
        Re = try
            stderror(fit)[2]
        catch
            0.
        end
        R, Re
    catch
        0., 0.
    end

    return R ± Re
end


function optimisewidth!(state)
    # make a range of dx values
    # set them, fit, and store the result
    # find the integration width giving minimum error in parameters
    # set the integration width to that value

    dx0 = state[:dx][]
    dxs = logrange(.02, 50, 21) .* dx0
    dxs = filter(x -> 0.01 <= x <= 5, dxs)

    # scoring function
    score(state) = try
        det(estimate_covar(state[:fit][]))
    catch e
        Inf
    end
    scores = map(dxs) do dx
        state[:dx][] = dx
        @debug dx, score(state)
        score(state)
    end

    best = argmin(scores)
    @debug "best dx: $(dxs[best])"
    state[:dx][] = dxs[best]
end