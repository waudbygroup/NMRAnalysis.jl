function linedata(spectra, plotscale)
    # produce a list of lists of Point2f containing spectrum data
    # x values are obtained from data(spec, F1Dim)
    # y values are obtained from data(spec) / scale
    
    map(enumerate(spectra)) do (idx, spec)
        x = data(spec, F1Dim)
        y = data(spec) * plotscale
        [Point2f(x[i], y[i] + idx - 1) for i in 1:length(x)]
    end
end

function linecolors(n)
    # using Colors

    # return [hue2rgb(HSV(i / n, 1.0, 1.0)) for i in 0:n-1]
    return :rainbow1
end