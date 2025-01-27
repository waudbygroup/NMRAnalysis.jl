function SpecData(nmrdata, x, y, z, σ, zlabels)
    zfit = Observable([0zi for zi in z])
    SpecData(nmrdata, x, y, z, σ, zlabels, zfit)
end