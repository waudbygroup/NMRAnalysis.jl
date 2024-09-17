# Model code: fitting amplitudes, keep everything else fixed

struct ModelAmplitudes <: Model end

# flags to control whether positions, amplitudes and linewidths are fixed across spectra
fixpeakpositions(::ModelAmplitudes) = true
fixamplitudes(::ModelAmplitudes) = false
fixlinewidths(::ModelAmplitudes) = true

function createpeak!(state, positionpoint, newlabel, ::ModelAmplitudes)
    # TODO set up initial parameters, dictionaries etc.
    Peak2D(positionpoint, newlabel)
end

function fitcluster!(state, cluster, ::ModelAmplitudes)
    @debug "fitting cluster (ModelAmplitudes)" cluster
    peaks = state[:peaks][] # for convenience
    specdata = state[:specdata]
    nspec = specdata[:nspec]
    Xradius = state[:Xradius][]
    Yradius = state[:Yradius][]

    # default parameters
    R2X0 = 30.
    R2Y0 = 30.
    minR2 = 3.
    maxR2 = 200.
    minamp = -1e6
    maxamp = 1e6

    # 1. get peak parameters
    npeaks = length(cluster)
    peakX = [peaks[peakindex].initial_position[1][1] for peakindex in cluster]
    peakY = [peaks[peakindex].initial_position[1][2] for peakindex in cluster]

    # 2. get spectrum region covering peaks in cluster
    Xmin = minimum(peakX) - Xradius
    Xmax = maximum(peakX) + Xradius
    Ymin = minimum(peakY) - Yradius
    Ymax = maximum(peakY) + Yradius
    @debug "ROI limits" Xmin Xmax Ymin Ymax

    # 3. extract ROI
    iXs = [Xmin .≤ specdata[:x][i] .≤ Xmax for i in 1:nspec]
    iYs = [Ymin .≤ specdata[:y][i] .≤ Ymax for i in 1:nspec]
    roispectra = [specdata[:spectra][i][Xmin .. Xmax, Ymin .. Ymax] for i in 1:nspec]
    Xroi = [specdata[:x][i][iXs[i]] for i in 1:nspec]
    Yroi = [specdata[:y][i][iYs[i]] for i in 1:nspec]
    Zroi = [specdata[:z][i][iXs[i], iYs[i]] for i in 1:nspec]

    # 4. get mask and mask data
    masks = [zeros(Bool, size(Zroi[i])) for i=1:nspec]
    for i = 1:nspec
        for j = 1:npeaks
            maskellipse!(masks[i], Xroi[i], Yroi[i], peakX[j], peakY[j], Xradius, Yradius)
        end
    end
    maskedZ = [Zroi[i][masks[i]] for i=1:nspec]
    # maskedZ = [Zroi[i] for i=1:nspec]

    # 5. initialise parameters
    # peakpars = [peakX ;; peakY ;; R2X ;; R2Y ;; amp...]
    peakpars = zeros(npeaks, 4+nspec)
    peakpars[:,1] .= peakX
    peakpars[:,2] .= peakY
    peakpars[:,3] .= R2X0
    peakpars[:,4] .= R2Y0
    for i = 1:nspec
        for j = 1:npeaks
            amp = roispectra[i][Near(peakX[j]), Near(peakY[j])]
            peakpars[j,4+i] = 1000 * amp
        end
    end

    minpars = 0peakpars
    minpars[:,1] = peakpars[:,1] .- Xradius
    minpars[:,2] = peakpars[:,2] .- Yradius
    minpars[:,3:4] .= minR2
    minpars[:,5:end] .= minamp

    maxpars = 0peakpars
    maxpars[:,1] = peakpars[:,1] .+ Xradius
    maxpars[:,2] = peakpars[:,2] .+ Yradius
    maxpars[:,3:4] .= maxR2
    maxpars[:,5:end] .= maxamp

    # 6. create pack/unpack closures
    # pack/unpack functions - fit everything, including chemical shifts
    pack(pars) = vec(pars)
    unpack(parvec) = reshape(parvec, :, 4+nspec)

    # pack initial parameters into vector
    p0 = pack(peakpars)
    pmin = pack(minpars)
    pmax = pack(maxpars)
    p0[p0 .< pmin] .= pmin[p0 .< pmin]
    p0[p0 .> pmax] .= pmax[p0 .> pmax]

    # 7. allocate space for simulated spectra
    Zsim = [zeros(size(Zroi[i])) for i=1:nspec]

    # 8. create simulation closure
    function sim(pars, getresid=false)
        # unpack parameter list
        pX = pars[:,1]
        pY = pars[:,2]
        pR2X = pars[:,3]
        pR2Y = pars[:,4]
        pamp = pars[:,5:end]

        for i=1:nspec
            Zsim[i] .= 0
            xaxis = dims(roispectra[i], F1Dim)
            yaxis = dims(roispectra[i], F2Dim)
            for j=1:npeaks
                simx = lineshape(xaxis, pX[j], pR2X[j])
                simy = lineshape(yaxis, pY[j], pR2Y[j])
                Zsim[i] += pamp[j,i] * simx .* simy'
            end
        end
        
        if getresid
            r = [vec(maskedZ[i] - Zsim[i][masks[i]]) for i=1:nspec]
            r = vcat(r...)
            return r
        else
            return Zsim
        end
    end

    # construct residuals function
    resid(p) = sim(unpack(p), true)

    # 9. run fit
    @debug "Running fit..."
    fit = LsqFit.lmfit(resid, p0, Float64[], #show_trace=true,
            autodiff=:finiteforward, #maxIter=50, 
            # x_tol=1e-3, g_tol=1e-6,
            lower=pmin, upper=pmax)
    pfit = unpack(coef(fit))
    @debug pfit

    # 10. update peak parameters
    for i=1:npeaks
        peakindex = cluster[i]
        peaks[peakindex].pars[:position] = Point2f(pfit[i, 1], pfit[i, 2])
        peaks[peakindex].pars[:R2X] = pfit[i, 3]
        peaks[peakindex].pars[:R2Y] = pfit[i, 4]
        peaks[peakindex].pars[:amplitude] = pfit[i, 5:end]
    end

    # 11. get parameter uncertainties
    # TODO

    # 12. update peak parameters for simulation
    for i=1:npeaks
        peakindex = cluster[i]
        peaks[peakindex].sim_parameters[:position] = MaybeVector(Point2f(pfit[i, 1], pfit[i, 2]))
        peaks[peakindex].sim_parameters[:R2X] = MaybeVector(pfit[i, 3])
        peaks[peakindex].sim_parameters[:R2Y] = MaybeVector(pfit[i, 4])
        peaks[peakindex].sim_parameters[:amplitude] = MaybeVector(pfit[i, 5:end])
    end

    # peaks = state[:peaks][]

    # for peakindex in cluster
    #     # copy 'results' direct from the initial position (TODO)
    #     newpos = deepcopy(peaks[peakindex].initial_position)
    #     ix = findnearest(state[:x][], newpos[1][1])
    #     iy = findnearest(state[:y][], newpos[1][2])
    #     R2X = 30.
    #     R2Y = 30.
    #     amp = state[:z][][ix, iy] * R2X * R2Y
        
    #     peaks[peakindex].sim_parameters[:position] = newpos
    #     peaks[peakindex].sim_parameters[:amplitude] = MaybeVector(amp)
    #     peaks[peakindex].sim_parameters[:R2X] = MaybeVector(R2X)
    #     peaks[peakindex].sim_parameters[:R2Y] = MaybeVector(R2Y)
    # end
end



