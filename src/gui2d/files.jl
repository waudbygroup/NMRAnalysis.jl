# format of peaks files:
# comments begin with #
# each line is a peak, with fields separated by whitespace
# fields are: label, xposition, yposition,

function loadpeaks!(expt)
    file = pick_file(; filterlist="peaks;old")
    file == "" && return

    @info "Loading peak file $file"
    isfitting = expt.isfitting[]
    if isfitting
        expt.isfitting[] = false
    end
    deleteallpeaks!(expt)
    readpeaklist!(expt, file)
    return expt.isfitting[] = isfitting
end

function saveresults!(expt)
    folder = pick_folder()
    folder == "" && return

    @info "Saving results to $folder"
    @async begin
        expt.state[][:mode][] = :fitting
        sleep(0.1) # allow time for mode change to be processed
    end
    @async begin # do fitting in a separate task
        sleep(0.2) # allow time for mode change to be processed

        writepeaklists!(expt, folder)

        # save post-fit parameters to separate file
        writefitresults!(expt, folder)

        # remove any files peak_XXX.pdf from output folder before saving new plots
        for file in readdir(folder)
            if occursin(r"^peak_.*\.pdf$", file)
                rm(joinpath(folder, file))
            end
        end
        save_peak_plots!(expt, folder)

        # save current contour plot figure
        save_contour_plot!(expt, folder)

        expt.state[][:mode][] = :normal
    end
end

"""
    read_peaklist!(expt, filepath::AbstractString)

Read a peak list file and add peaks to the experiment object. The file should contain
peaks in the format: 'label x_position y_position [optional_fields...]' with fields
separated by whitespace. Lines beginning with '#' are treated as comments.

# Arguments
- `expt`: Experiment object that implements `addpeak!`
- `filepath`: Path to the peak list file

# Returns
- Number of peaks successfully read

# Throws
- `ArgumentError`: If file format is invalid
- `IOError`: If file cannot be opened or read
"""
function readpeaklist!(expt, filepath::AbstractString)
    peak_count = 0

    open(filepath) do f
        for (line_number, line) in enumerate(eachline(f))
            # Skip empty lines and comments
            isempty(strip(line)) && continue
            startswith(strip(line), '#') && continue

            try
                fields = split(strip(line))
                if length(fields) < 3
                    throw(ArgumentError("Line $line_number: Insufficient fields"))
                end

                label = string(fields[1])  # Ensure string type
                x = parse(Float64, fields[2])
                y = parse(Float64, fields[3])

                addpeak!(expt, Point2f(x, y), label)
                peak_count += 1

            catch e
                if e isa ArgumentError
                    @warn "Skipping line $line_number: $(e.msg)"
                else
                    rethrow()
                end
            end
        end
    end

    @debug "Added $peak_count peaks from $filepath"
    return peak_count
end

"""
   write_peaklists!(expt::FixedPeakExperiment, folder) -> Vector{String}

Write peak lists for a FixedPeakExperiment.
"""
function writepeaklists!(expt::FixedPeakExperiment, folder)
    filepath_fit = joinpath(folder, "fit.peaks")
    writepeaklist!(filepath_fit, expt.peaks, 1, :value, expt)

    filepath_initial = joinpath(folder, "initial.peaks")
    return writepeaklist!(filepath_initial, expt.peaks, 1, :initialvalue, expt)
end

"""
   write_peaklists!(expt::MovingPeakExperiment, folder) -> Vector{String}

Write peak lists for a variable peak experiment.
"""
function writepeaklists!(expt::MovingPeakExperiment, folder)
    for slice in 1:nslices(expt)
        filepath_fit = joinpath(folder, "fit.$slice.peaks")
        writepeaklist!(filepath_fit, expt.peaks, slice, :value, expt)

        filepath_initial = joinpath(folder, "initial.$slice.peaks")
        writepeaklist!(filepath_initial, expt.peaks, slice, :initialvalue, expt)
    end
end

"""
   writepeaklist!(filepath, peaks, slice, value_type, expt) -> String

Write a single peak list file with specified format.
"""
function writepeaklist!(filepath, peaks, slice, value_type, expt)
    backup_file(filepath)

    open(filepath, "w") do f
        writeheader!(f, value_type, expt)
        return writepeaks!(f, peaks, slice, value_type, expt)
    end

    return filepath
end

"""
    backup_file(filepath::AbstractString)

Create a backup of an existing file by appending '.old'.
"""
function backup_file(filepath::AbstractString)
    isfile(filepath) || return
    backup_path = filepath * ".old"
    @debug "Backing up $filepath to $backup_path"
    return mv(filepath, backup_path; force=true)
end

"""
    format_param(peak, param, slice, value_type) -> String

Format a parameter value, returning "NA" if not found.
"""
function format_param(peak, param, slice, value_type)
    haskey(peak.parameters, param) || return "NA"
    val = getproperty(peak.parameters[param], value_type)[][slice]
    return string(to_value(val)) # convert from Observables to plain values
end

"""
   write_header!(io, value_type, expt::FixedPeakExperiment)

Write header for fixed peak experiments, including all amplitudes in one file.
"""
function writeheader!(io, value_type, expt::FixedPeakExperiment)
    for line in split(experimentinfo(expt), '\n')
        println(io, "# ", line)
    end

    header = ["label", "x", "y", "R2x", "R2y"]
    append!(header, ["amp$i" for i in 1:nslices(expt)])

    if value_type == :value
        append!(header, ["x_err", "y_err", "R2x_err", "R2y_err"])
        append!(header, ["amp$(i)_err" for i in 1:nslices(expt)])
    end

    return println(io, "# ", join(header, "\t"))
end

"""
   write_header!(io, value_type, expt::MovingPeakExperiment)

Write header for moving peak experiments, with single amplitude per file.
"""
function writeheader!(io, value_type, expt::MovingPeakExperiment)
    for line in split(experimentinfo(expt), '\n')
        println(io, "# ", line)
    end

    header = ["label", "x", "y", "R2x", "R2y", "amp"]

    if value_type == :value
        append!(header, ["x_err", "y_err", "R2x_err", "R2y_err", "amp_err"])
    end

    return println(io, "# ", join(header, "\t"))
end

"""
   write_peaks!(io, peaks, slice, value_type, expt::FixedPeakExperiment)

Write peaks for fixed peak experiments, including all amplitudes in one file.
"""
function writepeaks!(io, peaks, slice, value_type, expt::FixedPeakExperiment)
    basic_params = [:x, :y, :R2x, :R2y]

    for peak in peaks[]
        values = [peak.label[]]

        # Basic parameters
        append!(values,
                [format_param(peak, param, 1, value_type) for param in basic_params])

        # All amplitudes
        append!(values, [format_param(peak, :amp, i, value_type) for i in 1:nslices(expt)])

        # Uncertainties if fit result
        if value_type == :value
            append!(values,
                    [format_param(peak, param, 1, :uncertainty) for param in basic_params])
            append!(values,
                    [format_param(peak, :amp, i, :uncertainty) for i in 1:nslices(expt)])
        end

        println(io, join(values, "\t"))
    end
end

"""
   write_peaks!(io, peaks, slice, value_type, expt::MovingPeakExperiment)

Write peaks for moving peak experiments, with single amplitude per file.
"""
function writepeaks!(io, peaks, slice, value_type, expt::MovingPeakExperiment)
    basic_params = [:x, :y, :R2x, :R2y, :amp]

    for peak in peaks[]
        values = [peak.label[]]

        # Basic parameters and amplitude for this slice
        append!(values,
                [format_param(peak, param, slice, value_type) for param in basic_params])

        # Uncertainties (if a fit result)
        if value_type == :value
            append!(values,
                    [format_param(peak, param, slice, :uncertainty)
                     for param in basic_params])
        end

        println(io, join(values, "\t"))
    end
end

"""
    writepeaklist!(filepath, peaks, slice, value_type) -> String

Write a single peak list file with specified format.
"""
function writepeaklist!(filepath, peaks, slice, value_type)
    backup_file(filepath)

    open(filepath, "w") do f
        write_header!(f, value_type)
        return write_peaks!(f, peaks, slice, value_type)
    end

    return filepath
end

"""
    write_header!(io, value_type)

Write the parameter header line to the file.
"""
function write_header!(io, value_type)
    param_names = [:x, :y, :R2x, :R2y, :amp]
    header = ["label";
              param_names;
              value_type == :value ? [Symbol("$(p)_err") for p in param_names] : []]
    return println(io, "# ", join(header, "\t"))
end

"""
    write_peaks!(io, peaks, slice, value_type)

Write all peaks for a given slice to the file.
"""
function write_peaks!(io, peaks, slice, value_type)
    param_names = [:x, :y, :R2x, :R2y, :amp]

    for (peak_idx, peak) in enumerate(peaks[])
        values = [peak.label[]]

        # Add parameter values
        append!(values,
                [format_param(peak, param, slice, value_type) for param in param_names])

        # Add uncertainties for fit results
        if value_type == :value
            append!(values,
                    [format_param(peak, param, slice, :uncertainty)
                     for param in param_names])
        end

        println(io, join(values, "\t"))
    end
end

"""
    write_fixed_peaklist!(expt, folder) -> Vector{String}

Write peak lists for a FixedPeakExperiment.
"""
function write_fixed_peaklist!(expt, folder)
    writepeaklist!(joinpath(folder, "fit.peaks"), expt.peaks, 1, :value)
    return writepeaklist!(joinpath(folder, "initial.peaks"), expt.peaks, 1, :initialvalue)
end

"""
    write_variable_peaklist!(expt, folder) -> Vector{String}

Write peak lists for a variable peak experiment.
"""
function write_variable_peaklist!(expt, folder)
    for slice in 1:nslices(expt)
        for (filename, value_type) in [("fit.$slice.peaks", :value),
                                       ("initial.$slice.peaks", :initialvalue)]
            filepath = joinpath(folder, filename)
            push!(created_files, writepeaklist!(filepath, expt.peaks, slice, value_type))
        end
    end
end

function writefitresults!(expt, folder)
    filepath = joinpath(folder, "fit-results.txt")
    backup_file(filepath)

    open(filepath, "w") do f
        for line in split(experimentinfo(expt), '\n')
            println(f, "# ", line)
        end

        # Write header
        header = ["label"]
        param_names = collect(keys(first(expt.peaks[]).postparameters))
        for param in values(first(expt.peaks[]).postparameters)
            param_label = replace(param.label, " " => "_")
            append!(header, ["$(param_label)_value", "$(param_label)_uncertainty"])
        end
        println(f, "# ", join(header, "\t"))

        # Write one line per peak with all parameters
        for peak in expt.peaks[]
            values = [peak.label[]]
            for param in param_names
                push!(values, string(peak.postparameters[param].value[][1]))
                push!(values, string(peak.postparameters[param].uncertainty[][1]))
            end
            println(f, join(values, "\t"))
        end
    end
end

function save_contour_plot!(expt, folder)
    CairoMakie.activate!()

    state = expt.state[]
    rect = state[:gui][][:axcontour].finallimits[]
    x0, y0 = rect.origin
    dx, dy = rect.widths
    lims = ((x0, x0 + dx), (y0, y0 + dy))

    for i in 1:nslices(expt)
        f = Figure()
        ax = Axis(f[1, 1];
                  xlabel="$(label(expt.specdata.nmrdata[1],F1Dim)) chemical shift (ppm)",
                  ylabel="$(label(expt.specdata.nmrdata[1],F2Dim)) chemical shift (ppm)",
                  xreversed=true, yreversed=true, limits=lims,
                  title=slicelabel(expt, i))
        heatmap!(ax,
                 expt.specdata.x[i],
                 expt.specdata.y[i],
                 expt.specdata.mask[][i];
                 colormap=[:white, :lightgoldenrod1],
                 colorrange=(0, 1))
        contour!(ax,
                 expt.specdata.x[i],
                 expt.specdata.y[i],
                 expt.specdata.z[i];
                 levels=state[:gui][][:contourlevels],
                 color=bicolours(:grey50, :lightblue))
        contour!(ax,
                 expt.specdata.x[i],
                 expt.specdata.y[i],
                 expt.specdata.zfit[][i];
                 levels=state[:gui][][:contourlevels],
                 color=bicolours(:orangered, :dodgerblue))
        scatter!(ax, state[:positions]; markersize=10, marker=:x, color=:black)
        text!(ax, state[:positions]; text=state[:labels],
              fontsize=14,
              # font=:bold,
              offset=(8, 0),
              align=(:left, :center),
              color=:black)

        save(joinpath(folder, "spectrum-fit-$i.pdf"), f)
    end

    return GLMakie.activate!()
end