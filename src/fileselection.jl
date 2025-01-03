last_folder = Ref{String}("")

function select_expts(directory_path=""; title_filters=String[])
    if isempty(directory_path)
        directory_path = last_folder[]
    end
    if isempty(directory_path)
        directory_path = pwd()
    end

    # Pick a folder from the given directory path
    folder = pick_folder(directory_path)
    
    # Return an empty array if no folder is selected
    if isempty(folder)
        return []
    end
    last_folder[] = folder
    
    # If the selected folder is an experiment, return it in an array
    if isexpt(folder)
        return [folder]
    else
        # Otherwise, pick experiment folders within the selected folder
        return pick_expt_folders(folder; title_filters=title_filters)
    end
end

function pick_expt_folders(directory_path; title_filters=String[])
    # Get all experiment directories and their titles
    numbered_dirs = filter(ispath, readdir(directory_path; join=true))
    
    # Build menu options, skipping folders without title files
    options = String[]
    folder_paths = String[]
    for dir in numbered_dirs
        folder_num = basename(dir)

        title_path = joinpath(dir, "pdata", "1", "title")
        title_content = readline(title_path) |> strip
        
        # Apply title filters if any are specified
        if isempty(title_filters) || any(f -> occursin(f, title_content), title_filters)
            push!(options, "$(folder_num): $(title_content)")
            push!(folder_paths, dir)
        end
    end
    
    # Handle case where no folders match filters
    if isempty(options)
        # fall back to unfiltered options
        return pick_expt_folders(directory_path)
    end
    
    # Sort both arrays based on folder numbers
    perm = sortperm(options, by=x -> parse(Int, split(x, ':')[1]))
    options = options[perm]
    folder_paths = folder_paths[perm]
    
    # Create and display the multi-select menu
    menu = MultiSelectMenu(options)
    choices = request("Select folders to analyse:", menu)
    
    # Convert Set to Vector for indexing and return the full paths
    return folder_paths[collect(choices)]
end


function isexpt(directory)
    # check for existence of pdata/1/title file
    title_path = joinpath(directory, "pdata", "1", "title")
    isfile(title_path)
end
