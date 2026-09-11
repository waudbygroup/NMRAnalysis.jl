# Results file and region round-trip, mirroring `gui2d/files.jl`: `results.csv` is both the
# machine-readable output and the format a saved region list is restored from, so a
# multi-region kinetics session is reproducible rather than only screenshot-able.
#
# Layout (as in 2D): `#` comment lines carrying the experiment metadata, then a header row,
# then one row per (region × group). Regions travel as ordinary `region`/`lo`/`hi` columns,
# the way 2D peaks travel as `label`/`x`/`y`, and it is those columns `readregions!` reads
# back.

"""
    experimentinfo(expt) -> String

Multi-line description of the experiment, written as the comment header of `results.csv`.
Generic across every experiment - the analysis name, where the data came from, how many
spectra and which variables are arrayed are all already known to the framework - so
experiments need not override it (GUI2D's equivalent is written out per experiment).
"""
function experimentinfo(expt::Experiment1D)
    ds = dataset(expt)
    io = IOBuffer()
    println(io, "Analysis type: $(windowtitle(expt))")
    isempty(ds.label) || println(io, "Source: $(ds.label)")
    println(io, "Number of spectra: $(nplanes(ds))")
    if !isempty(ds.planes.vars)
        println(io, "Arrayed variables: $(join(keys(first(ds.planes.vars)), ", "))")
    end
    return String(take!(io))
end

"""
    valueerr(params, key) -> (value, uncertainty)

The two CSV cells for one parameter, as strings. `"NA"` where the parameter is absent, and
for the uncertainty of a quantity that carries none (a solvent viscosity looked up from
temperature is a plain number, where a fitted rate is a `Measurement`) - matching GUI2D's
`format_param`/`format_post`.
"""
function valueerr(params, key::Symbol)
    haskey(params, key) || return ("NA", "NA")
    v = params[key]
    v isa Measurement || return (string(v), "NA")
    return (string(Measurements.value(v)), string(Measurements.uncertainty(v)))
end

"""
    resultstable(expt, results, regions) -> (header, rows)

Column names and rows for `results.csv`: one row per (region × group), carrying the
region's bounds, the fit's own parameters and the quantities derived from them, with the
experiment's [`primaryparam`](@ref) first among the derived columns.

Columns are the union of the keys present across all results, so a table stays rectangular
when one region fitted and another did not; missing cells read `NA`.
"""
function resultstable(expt::Experiment1D, results, regs)
    bounds = Dict(r.label => r for r in regs)
    fitkeys = unique(Iterators.flatten(keys(r.parameters) for r in results))
    derived = unique(Iterators.flatten(keys(r.postparameters) for r in results))
    primary = primaryparam(expt)
    primary in derived && (derived = [primary; filter(!=(primary), derived)])

    header = ["region", "lo", "hi", "group"]
    for k in Iterators.flatten((fitkeys, derived))
        append!(header, [string(k), "$(k)_err"])
    end

    rows = Vector{String}[]
    for r in results
        reg = get(bounds, r.region, nothing)
        row = [r.region,
               isnothing(reg) ? "NA" : string(reg.lo),
               isnothing(reg) ? "NA" : string(reg.hi),
               groupname(r.group)]
        for k in fitkeys
            append!(row, collect(valueerr(r.parameters, k)))
        end
        for k in derived
            append!(row, collect(valueerr(r.postparameters, k)))
        end
        push!(rows, row)
    end
    return header, rows
end

"""
    seriescolumn(result) -> String

Column name for one series in `intensities.csv`: the region label, suffixed with the group
where an experiment has more than one series per region (`amide_trosy`).
"""
function seriescolumn(r::RegionResult)
    return isempty(r.group) ? r.region : "$(r.region)_$(groupname(r.group))"
end

"""
    intensitiestable(expt, results) -> (header, rows)

Column names and rows for `intensities.csv`: the reduced quantity of every series against
the evolution parameter. One row per value of the fit axis, one pair of columns (value and
uncertainty) per series, i.e. times down the side and regions across the top - which is
the shape a kinetics experiment wants, and readable straight into a spreadsheet or
`DataFrame` for plotting.

Written for every experiment, not just the unfitted ones: the integrals behind a fit are
worth having whether or not one was performed. For a `NoFitting` experiment it is the
*only* place the deliverable appears, `results.csv` carrying fitted parameters of which
such an experiment has none.

Series need not share a fit axis (separate kinetics runs generally do not), so the rows
are the sorted union of every series' own values and a series missing one reads `NA`
there.
"""
function intensitiestable(expt::Experiment1D, results)
    header = [string(fitaxis(expt))]
    isempty(results) && return (header, Vector{String}[])

    xs = sort(unique(Float64[x for r in results for x in r.x]))
    for r in results
        name = seriescolumn(r)
        append!(header, [name, "$(name)_err"])
    end

    rows = Vector{String}[]
    for x in xs
        row = [string(x)]
        for r in results
            i = findfirst(==(x), r.x)
            if isnothing(i)
                append!(row, ["NA", "NA"])
            else
                append!(row, [string(Measurements.value(r.y[i])),
                              string(Measurements.uncertainty(r.y[i]))])
            end
        end
        push!(rows, row)
    end
    return header, rows
end

"""
    writetable(filepath, expt, header, rows; extra=String[]) -> String

Write one CSV: the experiment description as `#` comment lines (plus any `extra` comment
lines), then `header`, then `rows`. Any existing file is backed up first. Returns the path.
"""
function writetable(filepath::AbstractString, expt::Experiment1D, header, rows;
                    extra=String[])
    backupfile(filepath)
    open(filepath, "w") do f
        for line in split(experimentinfo(expt), '\n')
            isempty(strip(line)) && continue
            println(f, "# ", line)
        end
        for line in extra
            println(f, "# ", line)
        end
        println(f, join(header, ","))
        for row in rows
            println(f, join(row, ","))
        end
    end
    return filepath
end

"""
    writeresults!(expt, state, folder) -> String

Write `results.csv` (fitted and derived parameters) and `intensities.csv` (the reduced
series they were fitted to) into `folder`, backing up any existing ones first. Returns the
path of `results.csv`.
"""
function writeresults!(expt::Experiment1D, state, folder::AbstractString)
    results, regs = state[:result][], state[:regions][]
    # The noise position is a user choice the region bounds do not capture, so it is
    # recorded for `readregions!` to restore.
    noise = ["Noise position / ppm: $(round(state[:noisec][]; digits=4))"]
    filepath = writetable(joinpath(folder, "results.csv"), expt,
                          resultstable(expt, results, regs)...; extra=noise)
    writetable(joinpath(folder, "intensities.csv"), expt,
               intensitiestable(expt, results)...; extra=noise)
    return filepath
end

"""
    readregions!(state, filepath) -> Int

Restore a saved region list (and the noise position) from a `results.csv`, replacing
whatever is currently set, and return the number of regions read. Regions are taken from
the `region`/`lo`/`hi` columns, de-duplicated by label and kept in file order - a grouped
experiment writes one row per (region, group), so the same region appears several times.

Nothing else is read back: the parameter columns are outputs, recomputed from the restored
regions the moment they are set.
"""
function readregions!(state, filepath::AbstractString)
    isfile(filepath) || throw(ArgumentError("no such results file: $filepath"))
    regs = Region[]
    seen = Set{String}()
    colmap = nothing
    for line in eachline(filepath)
        sline = strip(line)
        isempty(sline) && continue
        if startswith(sline, '#')
            m = match(r"Noise position / ppm:\s*(\S+)", sline)
            isnothing(m) || (state[:noisec][] = parse(Float64, m.captures[1]))
            continue
        end
        fields = strip.(split(sline, ','))
        if isnothing(colmap)
            colmap = Dict(lowercase(name) => i for (i, name) in enumerate(fields))
            all(haskey(colmap, c) for c in ("region", "lo", "hi")) ||
                throw(ArgumentError("$filepath has no region/lo/hi columns"))
            continue
        end
        label = String(fields[colmap["region"]])
        label in seen && continue
        lo, hi = fields[colmap["lo"]], fields[colmap["hi"]]
        (lo == "NA" || hi == "NA") && continue
        push!(seen, label)
        push!(regs, Region(label, parse(Float64, lo), parse(Float64, hi)))
    end
    isempty(regs) && throw(ArgumentError("$filepath contains no regions"))
    state[:regions][] = regs
    state[:active][] = 1
    return length(regs)
end

"""Rename an existing file to `<name>.bak` so a save never silently destroys the last one."""
function backupfile(filepath::AbstractString)
    isfile(filepath) || return nothing
    return mv(filepath, filepath * ".bak"; force=true)
end
