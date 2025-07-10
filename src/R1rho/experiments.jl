function processexperiments(experimentfiles, minνSL=250.0)
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

        offres = occursin("offres", expt[:pulseprogram])
        if offres
            ΩSL, νSL, TSL, spec = process_offres_experiment(expt)
            n = ones(length(spec))
        else
            ΩSL, νSL, TSL, spec = process_onres_experiment(expt)
            n = -1 * ones(length(spec))
        end
        append!(ΩSLs, ΩSL)
        append!(νSLs, νSL)
        append!(TSLs, TSL)
        append!(spectra, spec)
        append!(exptnumbers, n)
    end

    # normalise by maximum intensity
    mx = maximum(map(maximum, spectra))
    spectra ./= mx

    # remove low νSL from TSL, ΩSL, νSL and spectra lists
    idx = findall(νSLs .> minνSL)
    ΩSLs = ΩSLs[idx]
    νSLs = νSLs[idx]
    TSLs = TSLs[idx]
    spectra = spectra[idx]
    exptnumbers = exptnumbers[idx]

    return R1RhoDataset(exptnumbers, ΩSLs, νSLs, TSLs, spectra)
end

function process_offres_experiment(expt)
    # 2. Get list of spinlock offsets
    fqlist = acqus(expt, :fq1list)
    ΩSL = getoffset(fqlist, dims(expt, F1Dim)) # in Hz

    # 3. Get list of relaxation times
    TSL = acqus(expt, :vplist)

    # 4. Get spinlock power
    spinlock_power_W = acqus(expt, :plw, 25)
    νSL = convert_W_to_Hz(spinlock_power_W, expt)

    nΩ = length(ΩSL)
    nT = length(relaxation_times)

    νSL = vec([νSL for i in 1:nT, j in 1:nΩ])
    TSL = vec([TSL[i] for i in 1:nT, j in 1:nΩ])
    ΩSL = vec([ΩSL[j] for i in 1:nT, j in 1:nΩ])
    spec = vec([expt[:, j, i] for i in 1:nT, j in 1:nΩ])

    return ΩSL, νSL, TSL, spec
end

function process_onres_experiment(expt)
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

    return ΩSL, νSL, TSL, spec
end