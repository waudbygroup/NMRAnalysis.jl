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

function fitexp(t, I, p)
    I0 = p[1]
    R0 = p[3] / 2
    p0 = [I0, R0]

    model(t, p) = @. p[1] * exp(-t * p[2])

    fit = LsqFit.curve_fit(model, t, I, p0)
    R = coef(fit)[2]
    Re = try
        stderror(fit)[2]
    catch
        0.
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
    score(state) = det(estimate_covar(state[:fit][]))
    scores = map(dxs) do dx
        @info dx
        state[:dx][] = dx
        @show score(state)
    end

    best = argmin(scores)
    @info "best dx: $(dxs[best])"
    state[:dx][] = dxs[best]
end