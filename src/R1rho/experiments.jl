function processexperiments(experimentfiles)
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
        expt /= scale(expt)

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
    mx = map(maximum, spectra) |> maximum
    spectra ./= mx

    R1RhoDataset(exptnumbers, ΩSLs, νSLs, TSLs, spectra)
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

    νSL  = [νSL           for i in 1:nT, j in 1:nΩ] |> vec
    TSL  = [TSL[i]        for i in 1:nT, j in 1:nΩ] |> vec
    ΩSL  = [ΩSL[j]        for i in 1:nT, j in 1:nΩ] |> vec
    spec = [expt[:, j, i] for i in 1:nT, j in 1:nΩ] |> vec

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

    νSL  = [νSL[j]        for i in 1:nT, j in 1:nν] |> vec
    TSL  = [TSL[i]        for i in 1:nT, j in 1:nν] |> vec
    ΩSL  = [ΩSL           for i in 1:nT, j in 1:nν] |> vec
    spec = [expt[:, j, i] for i in 1:nT, j in 1:nν] |> vec

    return ΩSL, νSL, TSL, spec
end