function integrate(spectra::Vector, x1, x2)
    a,b = minmax(x1, x2)
    [integrate(s, a, b) for s in spectra]
end

function integrate(spectrum::NMRData, x1, x2)
    sum(spectrum[x1 .. x2])
end