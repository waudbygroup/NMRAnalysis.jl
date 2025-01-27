function SpecData(nmrdata, x, y, z, σ, zlabels)
    zfit = Observable([zeros(size(zi)) for zi in z])
    mask = Observable([falses(size(zi)) for zi in z])
    SpecData(nmrdata, x, y, z, σ, zlabels, zfit, mask)
end

