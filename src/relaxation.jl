function relaxation()
    println("Current directory: $(pwd())")
    println()

    print("Enter path to relaxation experiment (i.e. Bruker experiment folder): ")
    filename = readline()

    return relaxation(filename)
end

function relaxation(filename::String)
    ispath(filename) || throw(SystemError("No such file or directory"))

    print("Enter 'IR' for inversion-recovery otherwise press ENTER: ")
    input = readline()

    ir = lowercase(input) == "ir"
    return relaxation(loadnmr(filename); ir=ir)
end

function relaxation(spec::NMRData{T,2}; ir=false) where {T}
    # check for vdlist or vclist. If vclist, get conversion factor from user
    if haskey(acqus(spec), :vdlist) && acqus(spec, :vdlist) != ""
        tau = acqus(spec, :vdlist)
        println("Found vdlist in acqus. Using as relaxation delays.")
    elseif haskey(acqus(spec), :vclist)
        print("Found vclist in acqus. Please enter milliseconds per loop: ")
        vclist = acqus(spec, :vclist)
        conversion_factor = tryparse(Float64, readline())
        if conversion_factor === nothing
            throw(ArgumentError("Invalid conversion factor"))
        end
        tau = vclist .* conversion_factor / 1000
    else
        throw(KeyError("No vdlist or vclist found in acqus"))
    end
    return relaxation(spec, tau; ir=ir)
end

function relaxation(spec::NMRData{T,2}, tau; ir=false) where {T}
    @info tau
    spec = deepcopy(spec) # work on copy of the data
    label!(spec, "Relaxation")

    # 1. pick integration region
    plt = plot(spec[:, 1]; grid=true)
    hline!(plt, [0]; c=:grey)
    display(plt)

    print("Defining integration region - please enter first chemical shift: ")
    ppm1 = readline()
    ppm1 = tryparse(Float64, ppm1)
    print("Defining integration region - please enter second chemical shift: ")
    ppm2 = readline()
    ppm2 = tryparse(Float64, ppm2)

    ppm1, ppm2 = minmax(ppm1, ppm2)

    # 3. pick noise region
    vspan!(plt, [ppm1, ppm2]; label="integration region", alpha=0.2)
    display(plt)

    print("Enter a chemical shift in the center of the noise region: ")
    noiseppm = readline()
    noiseppm = tryparse(Float64, noiseppm)

    # create integration region and noise selectors

    ppmrange = ppm2 - ppm1
    noise1 = noiseppm - 0.5ppmrange
    noise2 = noiseppm + 0.5ppmrange

    roi = ppm1 .. ppm2
    noiseroi = noise1 .. noise2

    # plot region
    p2 = plot(spec[:, 1]; linecolor=:black)
    hline!(p2, [0]; c=:grey, primary=false)
    plot!(p2, spec[roi, 1]; fill=(0, :dodgerblue), linecolor=:navy,
          label="integration region")
    plot!(p2, spec[noiseroi, 1]; fill=(0, :orange), linecolor=:red, label="noise region",
          legend=:topright)
    title!(p2, "Integration regions")
    display(p2)
    print("Displaying integration and noise regions. Press enter to continue.")
    readline()

    # 4. integrate regions
    noise = vec(data(sum(spec[noiseroi, :]; dims=F1Dim)))
    noise = std(noise)

    integrals = vec(data(sum(spec[roi, :]; dims=F1Dim)))

    # normalise by max value
    noise /= maximum(integrals)
    integrals /= maximum(integrals)

    # 5. fit region
    model(t, p) = ir ? p[1] * (1 .- p[3] * exp.(-t * p[2])) : p[1] * exp.(-t * p[2])
    p0 = ir ? [1.0, 2 / maximum(tau), 2.0] : [1.0, 2 / maximum(tau)]

    fit = curve_fit(model, tau, integrals, p0)
    pars = coef(fit) .± stderror(fit)
    I0 = coef(fit)[1]
    rate = pars[2]

    # 7. plot fitted region
    x = LinRange(0, maximum(tau) * 1.1, 100)
    yfit = model(x, coef(fit))

    p1 = scatter(tau, (integrals .± noise) / I0; label="observed",
                 frame=:box,
                 xlabel="Relaxation time / s",
                 ylabel="Integrated signal",
                 title="Relaxation rate = $rate s⁻¹",
                 grid=nothing)
    plot!(p1, x, yfit / I0;
          label="fit",
          z_order=:back)
    hline!(p1, [0]; color=:black, lw=0, primary=false)

    println()
    message = """Relaxation results

Current directory: $(pwd())
Experiment: $(spec[:filename])

Integration region: $ppm1 - $ppm2 ppm
Noise region: $noise1 - $noise2 ppm

Fitted relaxation rate: $rate s⁻¹
Fitted relaxation time: $(1/rate) s
"""
    if ir
        message *= "Inversion-recovery amplitude: $(pars[3])\n"
    end
    @info message

    display(p1)

    # 8. make final plot

    println()
    print("Enter a filename to save figure (press enter to skip): ")
    outputfilename = readline()
    if length(outputfilename) > 0
        plot!(p1; title="")
        savefig(p1, outputfilename)
        println("Figure saved to $outputfilename.")
    end

    return p1
end