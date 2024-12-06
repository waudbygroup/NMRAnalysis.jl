function linedata(spectra, plotscale)
    # produce a list of lists of Point2f containing spectrum data
    # x values are obtained from data(spec, F1Dim)
    # y values are obtained from data(spec) / scale
    
    dat = map(enumerate(spectra)) do (idx, spec)
        x = data(spec, F1Dim)
        y = data(spec) * plotscale
        [Point2f(x[i], y[i] + idx - 1) for i in 1:length(x)]
    end

    # if length of dat is < 96, pad it with empty lists
    if length(dat) < 96
        append!(dat, [Point2f[] for i in 1:96 - length(dat)])
    end

    return dat
end

function linecolors(n)
    # using Colors

    c = resample_cmap(:darkrainbow, n)
    # pad to length 96
    if length(c) < 96
        append!(c, [c[1] for i in 1:96 - length(c)])
    end
    return c
    # return :rainbow1
end