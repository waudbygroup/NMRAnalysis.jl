function loadpeaks!(expt)
    file = pick_file(filterlist="peaks")
    file == "" && return

    @info "Loading peak file $file (TODO)"
    # TODO
end

function saveresults!(expt)
    folder = pick_folder()
    folder == "" && return

    @info "Saving results to $folder (TODO)"
    # TODO
end