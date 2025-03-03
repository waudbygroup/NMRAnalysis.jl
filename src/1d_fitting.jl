function fit_diffusion(spec::NMRData{T,2}, selector::Selector, plotrange=:, normalize=true;
	δ=nothing, Δ=nothing, σ=nothing, coherence=SQ(H1), solvent=nothing, temperature=nothing) where T
	fit_diffusion(spec, [selector], plotrange, normalize,
		coherence=coherence, δ=δ, Δ=Δ, σ=σ, solvent=solvent, temperature=temperature)
end

function fit_diffusion(spec::NMRData{T,2}, selectors::Vector{Selector}, plotrange=:, normalize=true;
		δ=nothing, Δ=nothing, σ=nothing, coherence=SQ(H1), solvent=nothing, temperature=nothing) where T
	any(isa.(dims(spec),GradientDimension)) || ArgumentError("spectrum has no gradient dimension - use setgradientlist() first")

	if isnothing(δ)
		δ = acqus(dat, :p, 30) * 2e-6  # gradient pulse length
		@info "assuming δ = 2*p30 = $δ s"
	end

	if isnothing(Δ)
		Δ = acqus(dat, :d, 20)  # diffusion delay
		@info "assuming Δ = d20 = $Δ s"
	end

	if isnothing(σ)
		σ = 0.9
		@info "assuming shape factor σ = $σ (suitable for SMSQ)"
	end

	if !isnothing(solvent)
		if isnothing(temperature)
			temperature = acqus(dat, :te)
			@info "reported temperature T = $temperature K"
		end
		η = viscosity(solvent, temperature)
	else
		η = nothing
	end
	
	γ = γeff(coherence)
	g = data(dims(spec, G1Dim))

	integrals = map(selectors) do selector
		if selector isa ArraySelector
			vec(data(sum(spec[selector,:],dims=F1Dim)))
		else
			vec(data(spec[selector,:]))
		end
	end

	# fitting
	model(g, p) = p[1] * exp.(-(γ*δ*σ*g).^2 .* (Δ - δ/3) .* p[2] .* 1e-10)
    p0 = [1.0, 1.0] # rough guess of scaled diffusion coefficient
	fits = map(integrals) do y
		curvefit(model, g, y / maximum(y), p0)
	end
	pars = map(first, fits)
	# errs = map(x->x[2], fits)
	D = map(fits) do fit
		(fit[1][2] ± fit[2][2]) * 1e-10
	end

	# plot results
	x = LinRange(0, maximum(g)*1.1, 50)
    yfits = map(p -> model(x, p), pars)

    p1 = plot(frame=:box,
			xlabel="G / T m-1",
			ylabel="Integrated signal",
			title="",
			ylims=(0,Inf),
			widen=true,
			grid=nothing)
	for i=1:length(selectors)
		scatter!(p1, g, integrals[i], label=string(selectors[i]))
		plot!(p1, x, yfits[i], primary=false)
	end

    p2 = plot(spec[plotrange,1],linecolor=:black)
	for (i, selector) in enumate(selectors)
    	plot!(p2, spec[selector,1], fill=0, linecolor=i)
	end
	hline!(p2, [0], c=:grey)
    title!(p2, "")
    #xlims!(p2, (0,10))

    plt = plot(p1, p2, layout=(2,1), size=(600,600))

	if isnothing(η)
		return D, plt
	else
		kB = 1.38e-23
		rH = kB*temperature / (6π*η*0.001 * D) * 1e9 # in nm
		return D, plt, rH
	end
end

	
function viscosity(solvent, T)
	if solvent==:h2o
        A = 802.25336
        a = 3.4741e-3
        b = -1.7413e-5
        c = 2.7719e-8
        gamma = 1.53026
        T0 = 225.334
    elseif solvent==:d2o
        A = 885.60402
        a = 2.799e-3
        b = -1.6342e-5
        c = 2.9067e-8
        gamma = 1.55255
        T0 = 231.832
    else
        @error "solvent not recognised (should be :h2o or :d2o)"
    end

    DT = T - T0
    k = 1.38e-23
	
    return A * (DT + a*DT^2 + b*DT^3 + c*DT^4)^(-gamma)
end



