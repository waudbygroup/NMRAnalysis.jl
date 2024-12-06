function loadfolder(foldername)
    folders = filter(x -> is_integer_string(x), readdir(foldername))
    sorted_folders = sort(folders, by = x -> parse(Int, x))
    
    rawspectra = map(folder -> loadnmr(joinpath(foldername, folder)), sorted_folders)

    spectra = map(rawspectra) do spec
        spec / scale(spec)
    end
    mx = maximum(map(maximum, spectra))
    spectra = map(s -> s / 0.2mx, spectra)


    return spectra
end

is_integer_string(s) = tryparse(Int, s) !== nothing

function saveresults(state)
    filename = state["resultsfilename"][]
    foldername = state["foldername"]
    # save a csv file with columns: ID (A1..A12..B1....H12), standard integral, unknown integral, ratio
    # include header rows with the foldername, and integration regions
    
    # get the data
    ids = [string(c, i) for i in 1:12, c in 'A':'H']
    standardintegrals = state["standardintegrals"][]
    unknownintegrals = state["unknownintegrals"][]
    ratios = state["integralratios"][]

    # open the file
    open(joinpath(foldername, filename), "w") do io
        println(io, "Folder: $foldername")
        println(io, state["gui"]["integration_label_1"][])
        println(io, state["gui"]["integration_label_2"][])
        println(io, "ID,Standard,Unknown,Ratio")
        for (id, si, ui, r) in zip(ids, standardintegrals, unknownintegrals, ratios)
            println(io, "$id,$si,$ui,$r")
        end
    end
end