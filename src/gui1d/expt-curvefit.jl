struct IntensityExperiment <: Experiment
    specdata::SpecData
    regions::Observable{Vector{Region}}
    state::Dict{Symbol,Observable}

    x # x values for the experiment - this could be a simple vector, or a matrix of there are multiple independent variables
    model::FittingModel
    visualisation::VisualisationStrategy

    function IntensityExperiment(specdata, model=NoFitting(), xvalues=nothing, visualisation=SimpleVisualisation())
        if isnothing(xvalues)
            n = sum(map(length, specdata.y))
            xvalues = range(1.0, n, n)
        end
        expt = new(specdata,
                Observable(Region[]),
                Dict{Symbol, Observable}(),
                xvalues,
                model,
                visualisation)
        setupexptobservables!(expt)
        preparestate!(expt)
        expt
    end
end

visualisationtype(expt::IntensityExperiment) = expt.visualisation

# DEFINE EXPERIMENTS ###########################################################

function pseudo2d(inputfilename)
    nmrdata = loadnmr(inputfilename)

    x = [data(spec, F1Dim)]
    y = [LinRange(1, size(nmrdata, 2), size(nmrdata, 2))]
    z = [1]
    xlabel = label(spec, F1Dim)
    ylabel = ""
    zlabels = [""]
    dat = [data(spec) / spec[:noise]]
    σ = [1.0]
    specdata = SpecData(nmrdata, dat, x, y, z, σ, xlabel, ylabel, zlabels)
    
    expt = IntensityExperiment(
        specdata,
        NoFitting()
    )
    
    # gui!(expt)
    expt
end

function diffusion(inputfilename, gmin, gmax, ramptype; δ=nothing, Δ=nothing, σ=nothing, Gmax=0.55, coherence=SQ(H1))
    nmrdata = loadnmr(inputfilename)
    
    if isnothing(δ)
        δ = acqus(nmrdata, :p, 30) * 2 * 1e-6
        @info "Diffusion: setting δ from acquisition file (2 x p30)" δ
    end
    if isnothing(Δ)
        Δ = acqus(nmrdata, :d, 20)
        @info "Diffusion: setting Δ from acquisition file (d20)" Δ
    end
    if isnothing(σ)
        gpnam = acqus(spec, :gpnam, 6)  # try to identify shape factor
        σ = 1
        if length(gpnam) ≥ 4
            if gpnam[1:4] == "SMSQ"
                σ = 0.9
            elseif gpnam[1:4] == "SINE"
                σ = 0.6366
            end
        end
        @info "Diffusion: setting σ from acquisition file (gpnam6)" σ
    end
    γ = gyromagneticratio(coherence)
    # TODO different ramp types
    g = LinRange(g1, g2, size(nmrdata, 2)) .* Gmax

    x = [data(spec, F1Dim)]
    y = [g]
    z = [1] # only a single experiment
    xlabel = label(spec, F1Dim)
    ylabel = "Gradient strength / T m⁻¹"
    zlabels = [""]
    dat = [data(spec) / spec[:noise]]
    σ = [1.0]
    specdata = SpecData(nmrdata, dat, x, y, z, σ, xlabel, ylabel, zlabels)
    
    expt = IntensityExperiment(
        specdata,
        DiffusionModel(δ, Δ, σ, γ),
        g,
        ModelFitVisualisation()
    )
    
    # gui!(expt)
    expt
end

# TODO - will need to introduce post-fitting step to calculate tau_c
function tract(trosyfilename, antitrosyfilename)
    trosy = loadnmr(trosyfilename)
    antitrosy = loadnmr(antitrosyfilename)

    trosy = setrelaxtimes(trosy, acqus(trosy, :vdlist), "s")
    antitrosy = setrelaxtimes(antitrosy, acqus(antitrosy, :vdlist), "s")

    trosytau = data(trosy, 2)
    antitrosytau = data(antitrosy, 2)

    x = [data(trosy, F1Dim), data(antitrosy, F1Dim)]
    y = [trosytau, antitrosytau]
    z = [1, 2]
    xlabel = label(trosy, F1Dim)
    ylabel = "Relaxation delay / s"
    zlabels = ["TROSY", "Anti-TROSY"]
    σ = [trosy[:noise] / scale(trosy), antitrosy[:noise] / scale(antitrosy)]
    dat = [data(trosy) / scale(trosy), data(antitrosy) / scale(antitrosy)]
    dat ./= σ[1]
    σ ./= σ[1]
    specdata = SpecData([trosy, antitrosy], dat, x, y, z, σ, xlabel, ylabel, zlabels)

    # # calculate field-dependent quantities
    # B0 = 2π * 1e6 * acqus(trosy, :bf1) / γH
    # ωN = 2π * 1e6 * acqus(trosy, :bf3)

    # p = μ0 * γH * γN * ħ / (8π * sqrt(2) * rNH^3)
    # c = B0 * γN * ΔδN / (3 * sqrt(2))
    # f = p * c * (3cos(θ)^2 - 1)

    # AFTER FITTING
    # # 6. calculate τc
    # ΔR = antitrosyR2 - trosyR2
    # ηxy = ΔR / 2

    # # based on analytical solution in Mathematica
    # # Solve[f (4*4/10*tc + (3*4/10*tc/(1 + (tc*\[Omega]N)^2))) == \[Eta], tc]
    # x2 = 21952 * f^6 * ωN^6 - 3025 * f^4 * ηxy^2 * ωN^8 + 625 * f^2 * ηxy^4 * ωN^10
    # x = sqrt(x2)
    # y3 = 1800 * f^2 * ηxy * ωN^4 + 125 * ηxy^3 * ωN^6 + 24 * sqrt(3) * x
    # y = cbrt(y3)
    # τc = (5 * ηxy) / (24 * f) -
    #      (336 * f^2 * ωN^2 - 25 * ηxy^2 * ωN^4) / (24 * f * ωN^2 * y) + y / (24 * f * ωN^2)
    # τc = 1e9 * τc

    # TODO form x values for model combining x/y
    
    expt = IntensityExperiment(
        specdata,
        NoFitting()
    )
    
    # gui!(expt)
    expt
end


# EXPERIMENT INFORMATION ######################################################

function experimentinfo(expt::IntensityExperiment)
    info = [
        "Analysis type: Intensity",
        "Model: $(typeof(expt.model))",
        "Filename: $(expt.specdata.nmrdata[1][:filename])",
        "Number of regions: $(length(expt.regions[]))",
        "Experiment title: $(expt.specdata.nmrdata[1][:title])"
    ]
    
    # Add model-specific information
    append!(info, modelinfo(expt.model, expt.x))
    
    join(info, "\n")
end

function slicelabel(expt::IntensityExperiment, idx)
    if length(expt.specdata.zlabels) == 1
        "Slice $idx of $(nslices(expt))"
    else
        "$(expt.specdata.zlabels[idx]) ($idx of $(nslices(expt)))"
    end
end

function regioninfo(expt::IntensityExperiment, idx)
    if idx == 0
        return "No region selected"
    end
    
    region = expt.regions[][idx]
    
    # Common region information
    info = [
        "Region: $(region.label[])",
        "",
        "$(region.x1[]) - $(region.x2[]) ppm",
    ]
    
    # Add model-specific parameters
    append!(info, modelparameterinfo(region, expt.model))
    
    join(info, "\n")
end

# FITTING ######################################################################

fit(integrals, expt::IntensityExperiment) = fit(integrals, expt, expt.model)

function fit(_, ::IntensityExperiment, ::NoFitting)
    Dict{Symbol, Parameter}()
end

function fit(integrals, expt::IntensityExperiment, model::FittingModel)
    @debug "Fitting model"
    x = expt.x
    y = integrals

    # Initial parameter estimates
    pdic = estimate_parameters(x, y, model)
    p0 = [pdic[par] for par in model.param_names]
    
    # Fit the model
    fit = curve_fit(model.func, x, y, p0)
    pfit = coef(fit)
    perr = stderror(fit)

    pars = Dict{Symbol, Parameter}()
    for (i, name) in enumerate(model.param_names)
        pars[Symbol(name)] = Parameter(name, pfit[i], perr[i])
    end

    return pars
end