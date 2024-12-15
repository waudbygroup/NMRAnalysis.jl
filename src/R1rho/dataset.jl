struct R1RhoDataset
    exptnumbers
    ΩSLs
    νSLs
    TSLs
    spectra
end

ΩSL(dataset) = sort(unique(dataset.ΩSLs))
νSL(dataset) = sort(unique(dataset.νSLs))
TSL(dataset) = sort(unique(dataset.TSLs))
nΩSL(dataset) = length(ΩSL(dataset))
nνSL(dataset) = length(νSL(dataset))
nTSL(dataset) = length(TSL(dataset))

isonres(dataset) = ΩSL(dataset) == [0.0]

function integrate(dataset, x, dx)
    integrationrange = (x-0.5dx) .. (x+0.5dx)
    map(dataset.spectra) do spectrum
        sum(spectrum[integrationrange])
    end
end

function noise(dataset, x, dx)
    integrationrange = (x-0.5dx) .. (x+0.5dx)
    y = map(dataset.spectra) do spectrum
        sum(spectrum[integrationrange])
    end
    std(y)
end

function onresseries(dataset)
    [[i for i in 1:length(dataset.νSLs) if dataset.νSLs[i] == νSL && dataset.exptnumbers[i] < 0] for νSL in νSL(dataset)]
end