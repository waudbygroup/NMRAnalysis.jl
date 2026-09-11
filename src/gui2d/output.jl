# Output tables, written to the layout and column rules in
# `docs/src/advanced/conventions.md`:
#
#   out/
#     summary.txt        the record to read
#     results.csv        one row per peak: positions, linewidths and derived parameters
#     series.csv         the measurements, one row per peak per plane
#     global.csv         parameters fitted once across every peak (a titration Kd)
#     summary.pdf
#     peaks/<label>.pdf, peaks/<label>.csv
#     cluster_*.pdf
#
# The **entity** is the peak. Anything that varies plane by plane - the amplitude always,
# and for a moving-peak experiment the positions and linewidths too - belongs in
# `series.csv`, keyed by the plane's own coordinate rather than by an `amp[7]` column whose
# meaning lived only in a comment line.

# ---- coordinates --------------------------------------------------------------

"""
    seriescoordinates(expt) -> Vector{Pair{Symbol,Any}}

What distinguishes one plane from another, as `name => values`. A vector gives one entry
per plane; a scalar is a setting constant across the experiment (a saturation time, a CPMG
relaxation delay) and is repeated down the rows.

Named per quantity rather than reduced to one `x` column, so that a saved file says what
its planes actually were. This is the same hook Exchange1D defines for the same reason.
"""
seriescoordinates(e::IntensityExperiment) = [coordinatename(e.model) => e.x]
seriescoordinates(e::MovingExperiment) = [coordinatename(e.model) => e.x]
function seriescoordinates(e::CESTExperiment)
    return [:offset => e.frequencies, :B1 => e.B1, :Tsat => e.Tsat]
end
seriescoordinates(e::CPMGExperiment) = [:nu_CPMG => e.vCPMG, :Trelax => e.Trelax]
seriescoordinates(e::HetNOEExperiment) = [:saturated => collect(e.saturation)]
seriescoordinates(e::CCRExperiment) = [:buildup => collect(e.isbuildup), :Trelax => e.T]

"""
    coordinatename(model) -> Symbol

What the per-plane `x` of an experiment fitted with `model` actually is. The generic
fallback is the uninformative `:x`, which is the honest answer for a model fitted against
whatever the caller supplied (`modelfit2d`).
"""
coordinatename(::FittingModel) = :x
coordinatename(::NoFitting) = :plane
coordinatename(::ExponentialModel) = :time
coordinatename(::RecoveryModel) = :time
coordinatename(::MethylCCRModel) = :time
coordinatename(::TitrationModel) = :concentration

const COORDINATE_UNITS = Dict(:time => "s", :offset => "ppm", :B1 => "Hz", :Tsat => "s",
                              :nu_CPMG => "Hz", :Trelax => "s")

coordinateunit(name::Symbol) = get(COORDINATE_UNITS, name, "")

# ---- units --------------------------------------------------------------------

# ASCII units for the fitted and derived parameters, as they appear in column headers. Only
# what can be stated with confidence: a parameter missing here gets no unit rather than a
# guessed one. `Kd` is the notable blank - a titration's concentrations come from the
# sample metadata, so its units are whatever that metadata used.
const PARAM_UNITS = Dict(:x => "ppm", :y => "ppm", :R2x => "s-1", :R2y => "s-1",
                         :R => "s-1", :R1 => "s-1", :R2 => "s-1", :R20 => "s-1",
                         :PRE => "s-1", :eta => "s-1", :S2tc => "ns",
                         :CSP => "ppm", :dX => "ppm", :dY => "ppm",
                         :Xfree => "ppm", :Xbound => "ppm",
                         :Yfree => "ppm", :Ybound => "ppm")

paramunit(name::Symbol) = get(PARAM_UNITS, name, "")

"""
    globalparams(expt) -> Vector{Symbol}

Derived parameters that are fitted once across every peak rather than per peak, and so
belong in `global.csv`. A titration's `Kd` is the case: it is fitted globally and then
copied onto every peak, so writing it in `results.csv` would repeat one number down a
column.
"""
globalparams(::Experiment) = Symbol[]
globalparams(e::MovingExperiment) = e.model isa TitrationModel ? [:Kd] : Symbol[]

# ---- provenance ---------------------------------------------------------------

"""Comment lines heading every CSV: the experiment description plus the fitting radii,
which `readpeaklist!` reads back."""
function resultcomments(expt)
    lines = ["NMRAnalysis.jl $(pkgversion(GUI2D))"]
    append!(lines, split(experimentinfo(expt), '\n'))
    push!(lines, "X radius / ppm: $(round(expt.xradius[]; digits=4))")
    push!(lines, "Y radius / ppm: $(round(expt.yradius[]; digits=4))")
    return lines
end

"""Where plane `i` came from. One file per plane for a multi-file experiment; the same
pseudo-3D dataset for every plane otherwise."""
function planesource(expt, i)
    nmrdata = expt.specdata.nmrdata
    spec = nmrdata[min(i, length(nmrdata))]
    return string(something(spec[:filename], ""))
end

# ---- results.csv --------------------------------------------------------------

"The lineshape parameters of a peak: where it is, and how broad it is in each dimension."
const POSITION_PARAMS = (:x, :y, :R2x, :R2y)


"""
    derivedkeys(expt) -> Vector{Symbol}

Derived parameter names for `results.csv`, the experiment's [`primaryparam`](@ref) first
and anything global left out.
"""
function derivedkeys(expt)
    peaks = expt.peaks[]
    isempty(peaks) && return Symbol[]
    keys_ = [k for k in keys(first(peaks).postparameters) if !(k in globalparams(expt))]
    primary = primaryparam(expt)
    i = findfirst(==(primary), keys_)
    isnothing(i) || pushfirst!(keys_, popat!(keys_, i))
    return keys_
end

"""
    resultstable(expt) -> (header, rows)

Column names and rows for `results.csv`: one row per peak, carrying its identity, its
position and linewidths where those are properties of the peak rather than of each plane,
and the parameters derived from the fit.

A moving-peak experiment's positions vary plane by plane and so are in `series.csv`
instead; only a fixed-peak experiment has a single position to report here.
"""
function resultstable(expt)
    peaks = sortedpeaks(expt)
    derived = derivedkeys(expt)
    fixed = hasfixedpositions(expt)

    header = ["label", "resnum", "resname", "atom"]
    if fixed
        for p in POSITION_PARAMS
            append!(header, collect(csvcolumns(p, paramunit(p))))
        end
    end
    for k in derived
        append!(header, collect(csvcolumns(k, paramunit(k))))
    end

    rows = Vector{String}[]
    for peak in peaks
        lbl = parse_label(peak.label[])
        row = [peak.label[], string(lbl.resnum),
               lbl.onelettercode == '?' ? "" : string(lbl.onelettercode), lbl.atom]
        if fixed
            for p in POSITION_PARAMS
                push!(row, format_param(peak, p, 1, :value))
                push!(row, format_param(peak, p, 1, :uncertainty))
            end
        end
        for k in derived
            push!(row, format_post(peak, k, :value))
            push!(row, format_post(peak, k, :uncertainty))
        end
        push!(rows, row)
    end
    return header, rows
end

# ---- series.csv ---------------------------------------------------------------

"""
    seriestable(expt) -> (header, rows)

Column names and rows for `series.csv`: one row per peak per plane, with the plane's
coordinates and the quantities measured there.

Every experiment contributes the amplitude. A moving-peak experiment contributes its
positions and linewidths too, those varying plane by plane; for a fixed-peak experiment
they are single values and stay in `results.csv`.
"""
function seriestable(expt)
    n = nslices(expt)
    coords = seriescoordinates(expt)
    values = hasfixedpositions(expt) ? (:amp,) : (:amp, POSITION_PARAMS...)

    header = ["source", "label"]
    append!(header, [csvcolumn(name, coordinateunit(name)) for (name, _) in coords])
    for k in values
        append!(header, collect(csvcolumns(k, paramunit(k))))
    end

    rows = Vector{String}[]
    for peak in sortedpeaks(expt), i in 1:n
        row = [planesource(expt, i), peak.label[]]
        for (_, value) in coords
            push!(row, csvvalue(value isa AbstractVector ? value[i] : value))
        end
        for k in values
            push!(row, format_param(peak, k, i, :value))
            push!(row, format_param(peak, k, i, :uncertainty))
        end
        push!(rows, row)
    end
    return header, rows
end

# ---- global.csv ---------------------------------------------------------------

"""
    globaltable(expt) -> (header, rows)

Column names and rows for `global.csv`: the parameters [`globalparams`](@ref) names, read
from the first peak that carries them, since a global fit gives every peak the same value.
"""
function globaltable(expt)
    header = ["parameter", "value", "error", "unit"]
    rows = Vector{String}[]
    peaks = expt.peaks[]
    for k in globalparams(expt)
        i = findfirst(p -> haskey(p.postparameters, k), peaks)
        isnothing(i) && continue
        push!(rows, [string(k), format_post(peaks[i], k, :value),
                     format_post(peaks[i], k, :uncertainty), paramunit(k)])
    end
    return header, rows
end

# ---- writing ------------------------------------------------------------------

"""
    writeresults!(expt, folder) -> String

Write `results.csv`, `series.csv`, `global.csv` (where the experiment fits anything
globally) and one `peaks/<label>.csv` per peak into `folder`, and return the path of
`results.csv`. The per-peak files hold that peak's own rows of `series.csv`, so the data
behind each plot sits beside it under the same basename.
"""
function writeresults!(expt, folder)
    comments = resultcomments(expt)
    filepath = writetable(joinpath(folder, "results.csv"), comments, resultstable(expt)...)

    header, rows = seriestable(expt)
    writetable(joinpath(folder, "series.csv"), comments, header, rows)

    labelcol = findfirst(==("label"), header)
    for peak in expt.peaks[]
        label = peak.label[]
        writetable(joinpath(folder, "peaks", "$(safename(label)).csv"), comments, header,
                   filter(row -> row[labelcol] == label, rows))
    end

    gheader, grows = globaltable(expt)
    isempty(grows) || writetable(joinpath(folder, "global.csv"), comments, gheader, grows)
    return filepath
end

"""
    writesummary(filepath, expt) -> String

Write `summary.txt`: the experiment description, the fitting radii, and the headline
parameter for every peak, rounded for reading rather than written at full precision.
"""
function writesummary(filepath, expt)
    backupfile(filepath)
    peaks = sortedpeaks(expt)
    primary = primaryparam(expt)
    open(filepath, "w") do f
        for line in resultcomments(expt)
            isempty(strip(line)) && continue
            println(f, line)
        end
        println(f)
        gheader, grows = globaltable(expt)
        if !isempty(grows)
            println(f, "Global parameters:")
            for row in grows
                println(f, "  $(row[1]): $(row[2]) +/- $(row[3]) $(row[4])")
            end
            println(f)
        end
        unit = paramunit(primary)
        println(f, "$(primary)$(isempty(unit) ? "" : " / $unit") by peak:")
        for peak in peaks
            haskey(peak.postparameters, primary) || continue
            value = tryparse(Float64, format_post(peak, primary, :value))
            err = tryparse(Float64, format_post(peak, primary, :uncertainty))
            isnothing(value) && continue
            rounded = round(value; sigdigits=4)
            suffix = isnothing(err) ? "" : " +/- $(round(err; sigdigits=2))"
            println(f, "  $(rpad(peak.label[], 12)) $rounded$suffix")
        end
        return nothing
    end
    return filepath
end
