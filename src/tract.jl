function analyse_tract(trosy::String, antitrosy::String)
    analyse_tract(
        loadnmr(trosy),
        loadnmr(antitrosy)
    )
end

function analyse_tract(trosy::NMRData{T,2}, antitrosy::NMRData{T,2}) where T
    # numerical constants
    μ0 = 4π * 1e-7
    γH = gyromagneticratio(H1)
    γN = gyromagneticratio(N15)
    ħ = 6.626e-34 / 2π
    rNH = 1.02e-10
    ΔδN = 160e-6
    θ = 17 * π / 180

    # 1. get experiment parameters
    trosy = deepcopy(trosy) # work on copies of the data
    antitrosy = deepcopy(antitrosy)
    label!(trosy, "TROSY")
    label!(antitrosy, "Anti-TROSY")
    trosy = setrelaxtimes(trosy, acqus(trosy, :vdlist), "s")
    antitrosy = setrelaxtimes(antitrosy, acqus(antitrosy, :vdlist), "s")

    trosytau = data(trosy, 2)
    antitrosytau = data(antitrosy, 2)

    B0 = 2π * 1e6 * acqus(trosy, :bf1) / γH
    ωN = 2π * 1e6 * acqus(trosy, :bf3)

    p = μ0 * γH * γN * ħ / (8π * sqrt(2) * rNH^3)
    c = B0 * γN * ΔδN / (3*sqrt(2))
    f = p * c * (3cos(θ)^2 - 1)
    
    # 2. pick noise position
    plt = plot([
        trosy[:,1],
        antitrosy[:,1]
    ])
    display(plt)
    print("Enter a chemical shift in the noise region: ")
    noiseppm = readline()
    noiseppm = tryparse(Float64, noiseppm)
    # closeall()

    # 3. pick region
    plt = plot([
        trosy[:,1],
        antitrosy[:,1]
    ])
    vline!([noiseppm], label="noise")
    display(plt)
    print("Defining integration region - please enter first chemical shift: ")
    ppm1 = readline()
    ppm1 = tryparse(Float64, ppm1)
    print("Defining integration region - please enter second chemical shift: ")
    ppm2 = readline()
    ppm2 = tryparse(Float64, ppm2)
    # closeall()
    
    # create integration region and noise selectors
    ppm1, ppm2 = minmax(ppm1, ppm2)

    ppmrange = ppm2 - ppm1
    noise1 = noiseppm - 0.5ppmrange
    noise2 = noiseppm + 0.5ppmrange

    roi = ppm1 .. ppm2
    noiseroi = noise1 .. noise2

    # plot region
    p2 = plot(trosy[:,1],linecolor=:black)
    plot!(p2, trosy[roi,1], fill=(0,:dodgerblue), linecolor=:navy)
    plot!(p2, trosy[noiseroi,1], fill=(0,:orange), linecolor=:red)
	hline!(p2, [0], c=:grey)
    title!(p2, "Integration regions")
    display(p2)
    print("Displaying integration region. Press enter to continue.")
    readline()
    


    # 4. integrate regions
    trosynoise = vec(data(sum(trosy[noiseroi,:], dims=F1Dim)))
    antitrosynoise = vec(data(sum(antitrosy[noiseroi,:], dims=F1Dim)))
	trosynoise = std(trosynoise)
	antitrosynoise = std(antitrosynoise)
	
	trosyintegrals = vec(data(sum(trosy[roi,:],dims=F1Dim)))
    antitrosyintegrals = vec(data(sum(antitrosy[roi,:],dims=F1Dim)))
	
	# normalise by max value
	trosynoise /= maximum(trosyintegrals)
    antitrosynoise /= maximum(antitrosyintegrals)

    trosyintegrals /= maximum(trosyintegrals)
    antitrosyintegrals /= maximum(antitrosyintegrals)


    # 5. fit region
    model(t, p) = p[1] * exp.(-p[2] * t)
    p0 = [1.0, 5.0] # initial guesses for amplitude, R2
    trosyfit = curve_fit(model, trosytau, trosyintegrals, p0)
    antitrosyfit = curve_fit(model, antitrosytau, antitrosyintegrals, p0)
    trosypars = coef(trosyfit) .± stderror(trosyfit)
    antitrosypars = coef(antitrosyfit) .± stderror(antitrosyfit)
    trosyR2 = trosypars[2]
    antitrosyR2 = antitrosypars[2]

    # 6. calculate τc
    ΔR = antitrosyR2 - trosyR2
    ηxy = ΔR / 2

    # based on analytical solution in Mathematica
    # Solve[f (4*4/10*tc + (3*4/10*tc/(1 + (tc*\[Omega]N)^2))) == \[Eta], tc]
    x2 = 21952*f^6*ωN^6 - 3025*f^4*ηxy^2*ωN^8 + 625*f^2*ηxy^4*ωN^10
    x = sqrt(x2)
    y3 = 1800*f^2*ηxy*ωN^4 + 125*ηxy^3*ωN^6 + 24*sqrt(3)*x
    y = cbrt(y3)
    τc = (5*ηxy)/(24*f) - (336*f^2*ωN^2 - 25*ηxy^2*ωN^4)/(24*f*ωN^2*y) + y/(24*f*ωN^2)
    τc = 1e9 * τc
    
    # 7. plot fitted region
    maxt = max(maximum(trosytau), maximum(antitrosytau))
    t = LinRange(0, maxt*1.1, 100)
    yfit1 = model(t, coef(trosyfit))
    yfit2 = model(t, coef(antitrosyfit))

    p1 = scatter(trosytau, trosyintegrals .± trosynoise, label="TROSY",
			frame=:box,
			xlabel="Relaxation time / s",
			ylabel="Integrated signal",
			title="TRACT: τc = $τc ns",
			grid=nothing)
    scatter!(p1, antitrosytau, antitrosyintegrals .± antitrosynoise, label="Anti-TROSY")
    plot!(p1, t, yfit1, label="TROSY fit (R₂ = $trosyR2 s⁻¹)", c=1)
    plot!(p1, t, yfit2, label="Anti-TROSY fit (R₂ = $antitrosyR2 s⁻¹)", c=2)

    println("Fitted TROSY relaxation rate: $trosyR2 ns")
    println("Fitted anti-TROSY relaxation rate: $antitrosyR2 ns")
    println("Fitted τc: $τc ns")
    display(p1)

    # 6. adjust region?
    # 7. add new region?
    # 8. make final plot

    print("Enter a filename to save figure (press enter to skip): ")
    outputfilename = readline()
    if length(outputfilename) > 0
        savefig(p1, outputfilename)
        println("Figure saved to $outputfilename.")
    end
end