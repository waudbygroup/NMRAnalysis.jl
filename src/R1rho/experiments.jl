function processexperiments(experimentfiles; minvSL=250.0, maxvSL=1e6)
    # prepare empty list of floats for ΩSL, νSL and R1rho
    exptnumbers = Vector{Int64}()
    ΩSLs = Vector{Float64}()
    νSLs = Vector{Float64}()
    TSLs = Vector{Float64}()
    spectra = Vector{NMRData}()

    i = 0 # counter for experiments
    for experimentfile in experimentfiles
        i = i + 1
        expt = loadnmr(experimentfile)
        expt /= NMRTools.scale(expt)

        offres = if isnothing(NMRTools.annotations(expt))
            occursin("offres", expt[:pulseprogram])
        else
            "off-resonance" in NMRTools.annotations(expt, "features")
        end
        if offres
            ΩSL, νSL, TSL, spec = process_offres_experiment(expt)
            n = ones(length(spec))
        else
            ΩSL, νSL, TSL, spec = process_onres_experiment(expt, minvSL, maxvSL)
            n = -1 * ones(length(spec))
        end
        append!(ΩSLs, ΩSL)
        append!(νSLs, νSL)
        append!(TSLs, TSL)
        append!(spectra, spec)
        append!(exptnumbers, n)
    end

    # check that we have at least one experiment
    if isempty(spectra)
        @error "No measurements selected from the provided files."
        return nothing
    end

    # normalise by maximum intensity
    mx = maximum(map(maximum, spectra))
    spectra ./= mx

    return R1RhoDataset(exptnumbers, ΩSLs, νSLs, TSLs, spectra)
end

function process_offres_experiment(expt)
    # 2. Get list of spinlock offsets
    fqlist = acqus(expt, :fq1list)
    ΩSL = hz(fqlist, dims(expt, F1Dim)) # in Hz

    # 3. Get list of relaxation times
    TSL = acqus(expt, :vplist)

    # 4. Get spinlock power (using 90 degree p1@pl1 as reference)
    spinlock_power = acqus(expt, :plw, 25)
    νSL = hz(spinlock_power, acqus(expt, :pl, 1), acqus(expt, :p, 1), 90)

    nΩ = length(ΩSL)
    nT = length(relaxation_times)

    νSL = vec([νSL for i in 1:nT, j in 1:nΩ])
    TSL = vec([TSL[i] for i in 1:nT, j in 1:nΩ])
    ΩSL = vec([ΩSL[j] for i in 1:nT, j in 1:nΩ])
    spec = vec([expt[:, j, i] for i in 1:nT, j in 1:nΩ])

    return ΩSL, νSL, TSL, spec
end

function process_onres_experiment(expt, minvSL, maxvSL)
    # NB. we expect a 3D experiment, chemical shift * offset * relaxation time

    # 2. Get list of spinlock strengths
    νSL = convert_valist_to_Hz(expt)

    # 3. Get list of relaxation times
    TSL = acqus(expt, :vplist)

    # 4. Get spinlock offset
    ΩSL = 0

    nν = length(νSL)
    nT = length(TSL)

    νSL = vec([νSL[j] for i in 1:nT, j in 1:nν])
    TSL = vec([TSL[i] for i in 1:nT, j in 1:nν])
    ΩSL = vec([ΩSL for i in 1:nT, j in 1:nν])
    spec = vec([expt[:, j, i] for i in 1:nT, j in 1:nν])

    # 5. Filter out unwanted spinlock strengths
    low_vals = unique(round.(νSL[νSL .< minvSL]; digits=1))
    high_vals = unique(round.(νSL[νSL .> maxvSL]; digits=1))
    if !isempty(low_vals)
        @info "Filtering out spinlock strengths below $minvSL Hz ($(expt[:filename])): $(low_vals) Hz"
    end
    if !isempty(high_vals)
        @info "Filtering out spinlock strengths above $maxvSL Hz ($(expt[:filename])): $(high_vals) Hz"
    end
    idx = findall(minvSL .< νSL .< maxvSL)
    ΩSL = ΩSL[idx]
    νSL = νSL[idx]
    TSL = TSL[idx]
    spec = spec[idx]

    return ΩSL, νSL, TSL, spec
end