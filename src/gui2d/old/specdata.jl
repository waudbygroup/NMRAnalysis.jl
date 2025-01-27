function preparespecdata(spectra, xrange=:, yrange=:)
    N = length(spectra)
    specs = [spectra[i][xrange, yrange]/spectra[i][:noise] for i=1:N] # TODO normalisation
    X = [data(spec, F1Dim) for spec in specs]
    Y = [data(spec, F2Dim) for spec in specs]
    Z = [data(spec) for spec in specs]
    xlabels = [label(spectra[i], F1Dim) for i=1:N]
    ylabels = [label(spectra[i], F2Dim) for i=1:N]
    titles = [choptitle(label(spectra[i])) for i=1:N]

    Zfit = [0z for z in Z]
    # Zfit = [z+randn(size(z)...) for z in Z]

    Dict(
        :nspec => N,
        :spectra => specs,
        :x => X,
        :y => Y,
        :z => Z,
        :xlabels => xlabels,
        :ylabels => ylabels,
        :titles => titles,
        :zfit => Zfit,
    )
end


function choptitle(title, maxlength=30)
    if length(title) > maxlength
        title[1:maxlength] * "…"
    else
        title
    end
end
