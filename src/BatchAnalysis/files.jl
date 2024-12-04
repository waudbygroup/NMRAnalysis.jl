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