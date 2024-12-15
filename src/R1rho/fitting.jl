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