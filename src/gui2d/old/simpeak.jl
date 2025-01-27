function lineshape(ω, R, ωobs, w::ExponentialWindow)
    x = @. -R + 1im*(ωobs - ω) - π*w.lb
    T = w.tmax
    return @. real((1 - exp(T*x)) / x)
end

function lineshape(ω, R, ωobs, w::CosWindow)
    x = @. -R + 1im*(ωobs - ω)
    T = w.tmax
    Tx = T * x
    return @. real(T * (π*exp(Tx) + 2*Tx) / (π^2 + 4*Tx^2))
end

function lineshape(ω, R, ωobs, w::Cos²Window)
    x = @. -R + 1im*(ωobs - ω)
    Tx = w.tmax * x
    return @. real((π^2*(1-exp(Tx)) + 2*Tx^2) / ((π^2 + Tx^2) * x))
end

Base.Broadcast.broadcastable(w::WindowFunction) = Ref(w)
  

"Return a matrix of simulated intensities"
function simpeaks(pamp, pX, pY, plwX, plwY, Xroi, Yroi, spec)
    Zsim = zeros(length(Xroi), length(Yroi))

    # convert chemical shift axes to angular frequencies
    bfX = spec[1, :bf]
    bfY = spec[2, :bf]
    ωobsX = reshape(2π * Xroi * bfX, :, 1)
    ωobsY = reshape(2π * Yroi * bfY, 1, :)

    # get window functions
    windowX = spec[1, :window]
    windowY = spec[2, :window]

    for i=1:length(pX)
        # convert peak positions to angular frequencies
        ωX = 2π * pX[i] * bfX
        ωY = 2π * pY[i] * bfY

        lineshapeX = 1e4 * pamp[i] * lineshape(ωX, π * plwX[i], ωobsX, windowX)
        lineshapeY = lineshape(ωY, π * plwY[i], ωobsY, windowY)

        Zsim += lineshapeX * lineshapeY
    end

    return Zsim
end

function simpeaks(state)
    #(peaks, dX, dY, spec)

    npeaks = length(state.peaks[])
    peakX = [state.peaks[][i].position[1] for i=1:npeaks]
    peakY = [state.peaks[][i].position[2] for i=1:npeaks]
    amp = [state.peaks[][i].amp for i=1:npeaks]
    lwX = [state.peaks[][i].lwX for i=1:npeaks]
    lwY = [state.peaks[][i].lwY for i=1:npeaks]

    if length(amp[1])>1
        nT = length(amp[1])
        Z = zeros(size(state.specdata.Z))
        for i=1:nT
            amp = [state.peaks[][j].amp[i] for j=1:npeaks]
            Z[:,:,i] = simpeaks(amp, peakX, peakY, lwX, lwY, state.specdata.dX, state.specdata.dY, state.specdata.spectrum)
        end
        return Z
    else
        simpeaks(amp, peakX, peakY, lwX, lwY, state.specdata.dX, state.specdata.dY, state.specdata.spectrum)
    end
end