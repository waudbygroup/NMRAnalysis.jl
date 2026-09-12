"""
    Experiment1D

Abstract supertype for 1D analyses. A concrete experiment is a thin composition that
supplies a dataset, a list of regions, a series model, and the
fit-axis / grouping designation. The generic [`analyse`](@ref) pipeline does the rest;
experiments derive further quantities in [`postfit!`](@ref) / [`postfitglobal!`](@ref),
and may override `analyse` entirely for non-curve-fit shapes.

# Adding an experiment

One file per experiment, `expt-<name>.jl`, included at the foot of this file, laid out in
a fixed five-section order (see `expt-tract.jl` for the fullest example):

1. **entry point** — the exported `<name>1d(...)` function: load the spectra, pull the
   acquisition parameters, build the experiment, `run1d(expt; integration)`.
2. **type** — the `struct <Name>Experiment <: Experiment1D` and its keyword constructor.
3. **interface** — only the hooks that differ from the defaults below.
4. **science** — the series model, the derived-quantity maths, physics constants.
5. **presentation** — `windowtitle`, `resultlabels`, `spectruminfo`, and friends.

Interface (with defaults):
- `dataset(e)`        — the `Dataset1D` (default: `e.dataset`)
- `regions(e)`        — `Vector{Region}`, length ≥ 1 (default: `e.regions`)
- `integrate(region, e)` — region × planes → one quantity per plane (default: integration)
- `seriesmodel(e)`    — `SeriesModel` (default: `e.model`)
- `fitaxis(e)`        — `Symbol` naming the evolution variable (required)
- `groupcols(e)`      — `Tuple` of grouping variables (default `()`)
- `postfit!(rs, e)`   — derived quantities from one region's series (default: none)
- `postfitglobal!(results, e)` — derived quantities spanning regions (default: none)
- `primaryparam(e)`   — the headline quantity (default `:A`)
"""
abstract type Experiment1D end

# Field-assuming defaults, in the style of GUI2D's `nslices(expt)`. An experiment whose
# fields are named `dataset`/`regions`/`model` needs none of these three methods; one
# that computes them (or names them differently) overrides just the one it needs.
dataset(e::Experiment1D) = e.dataset
regions(e::Experiment1D) = e.regions
seriesmodel(e::Experiment1D) = e.model

groupcols(::Experiment1D) = ()

"""
    integrate(region, expt[, dataset]) -> Vector{Measurement{Float64}}

Stage 1 of an analysis (see `docs/src/advanced/pipeline.md`): reduce `region` to one
measured quantity per plane. The 1D counterpart of GUI2D's lineshape `fit!(cluster, expt)`,
and like it dispatched on the experiment, so one that measures its regions some other way
overrides this method. The default integrates.

`dataset` is separate from the experiment's own so the GUI can measure against
interactively-positioned regions and noise.
"""
integrate(region::Region, e::Experiment1D, ds::Dataset1D=dataset(e)) = integrate(region, ds)

"""
    RegionResult

The analysis of one series: one region × one grouping key. The 1D analogue of GUI2D's
`Peak`, and like it the *uniform* container every experiment fills, however different
their science.

# Fields
- `region`    : region label
- `group`     : `NamedTuple` of grouping values (empty if ungrouped)
- `x`         : evolution-parameter values (sorted)
- `y`         : reduced quantities (`Vector{Measurement}`)
- `planes`    : which plane of the dataset each point came from, in the same order, so the
  series file can name its source spectrum
- `parameters`     : the curve fit's own parameters (`A`, `R`, `ν`, …)
- `postparameters` : quantities derived from them by [`postfit!`](@ref) /
  [`postfitglobal!`](@ref) (`tauc`, `pulse90`, `rH`, …)
- `postscopes`     : what each derived quantity describes - see [`setpost!`](@ref)
- `model`     : the series model used
- `converged` : whether the fit converged
- `postfitted`: whether a derived result has been recorded

Both parameter dictionaries are `Symbol`-keyed, mirroring `peak.parameters` /
`peak.postparameters`, so generic consumers (the summary formatter, `primaryparam`) read
them by name without knowing the experiment. Values are untyped because not every derived
quantity carries an uncertainty - a solvent viscosity looked up from temperature is a
plain number, where a fitted rate is a `Measurement`.

Derived quantities are stored **in the unit named by `PARAM_UNITS`** (a 90° pulse in µs,
τc in ns), not in SI, so that one stored number serves both the summary and any later
tabular export.
"""
mutable struct RegionResult
    region::String
    group::NamedTuple
    x::Vector{Float64}
    y::Vector{Measurement{Float64}}
    planes::Vector{Int}
    parameters::OrderedDict{Symbol,Any}
    postparameters::OrderedDict{Symbol,Any}
    postscopes::OrderedDict{Symbol,Symbol}
    model::Any
    converged::Bool
    postfitted::Bool
end

# Convenience constructor for the common case of a freshly fitted series, which has no
# derived quantities yet.
function RegionResult(region, group, x, y, planes, parameters, model, converged)
    return RegionResult(region, group, x, y, planes, parameters,
                        OrderedDict{Symbol,Any}(), OrderedDict{Symbol,Symbol}(), model,
                        converged, false)
end

"""
    param(result, name) -> value

Value of parameter `name` (a `Symbol` or a `String`) from a `RegionResult`: the fit's own
parameters first, then the quantities derived from them by [`postfit!`](@ref) /
[`postfitglobal!`](@ref). One accessor for both, so `param(r, :R)` and `param(r, :tauc)`
read alike and a caller need not know which stage produced a quantity.
"""
function param(r::RegionResult, name::Symbol)
    haskey(r.parameters, name) && return r.parameters[name]
    haskey(r.postparameters, name) && return r.postparameters[name]
    throw(KeyError(name))
end
param(r::RegionResult, name::AbstractString) = param(r, Symbol(name))

"""
    setpost!(result, name, value; scope=:series)

Record a derived quantity on `result` and mark it post-fitted. The `postfit!` /
`postfitglobal!` counterpart of writing into `peak.postparameters` and setting
`peak.postfitted[]` in GUI2D.

`scope` says what the quantity actually describes, which decides where it is written (see
`docs/src/advanced/conventions.md`):

- `:series` (the default) - this series alone, e.g. nutation's 90° pulse length. Written on
  this result's own row of `results.csv`.
- `:region` - the whole region, across every group. TRACT's τc combines the TROSY and
  anti-TROSY rates, so it describes the region rather than either component; `postfit!`
  still has to record it on one of the two, and the scope is what stops that choice
  leaking into the output. Written as a row with the group keys left blank.
- `:global` - the whole analysis, e.g. a solvent viscosity, which is a property of the
  sample and the same for every region. Written to `global.csv`.
"""
function setpost!(r::RegionResult, name::Symbol, value; scope::Symbol=:series)
    scope in (:series, :region, :global) ||
        throw(ArgumentError("unknown parameter scope :$scope (expected :series, :region or :global)"))
    r.postparameters[name] = value
    r.postscopes[name] = scope
    r.postfitted = true
    return value
end

"""
    postscope(result, name) -> Symbol

The scope recorded for a derived quantity by [`setpost!`](@ref), defaulting to `:series`
for anything written into `postparameters` directly.
"""
postscope(r::RegionResult, name::Symbol) = get(r.postscopes, name, :series)

"""
    postfit!(results, expt)

Stage 2: derive the quantities reported for one region from its fitted series, recording
them with [`setpost!`](@ref) - nutation's 90° pulse length from ν, diffusion's rH from D,
TRACT's τc from the TROSY and anti-TROSY rates together. The default does nothing
(relaxation and kinetics derive nothing).

`results` holds every series of a single region, so anything combining conditions belongs
here rather than in [`postfitglobal!`](@ref). Mirrors GUI2D's `postfit!(peak, expt)`, which
likewise sees the whole peak.
"""
postfit!(::AbstractVector{RegionResult}, ::Experiment1D) = nothing

"""
    postfitglobal!(results, expt)

Stage 3: fit or derive quantities spanning every region of the analysis, and record them on
the relevant results. Runs after every `postfit!`. Mirrors GUI2D's `postfitglobal!(expt)`.
"""
postfitglobal!(::AbstractVector{RegionResult}, ::Experiment1D) = nothing

"""
    primaryparam(expt) -> Symbol

The experiment's headline quantity: the parameter a reader wants first, listed first
among the derived columns of any tabular export. Defaults to the amplitude `A` that every
`CurveFitModel` here fits. Mirrors GUI2D's `primaryparam(expt)`.
"""
primaryparam(::Experiment1D) = :A

"""
    seriesresults(e, [dataset, regions]; isfitting=true) -> Vector{RegionResult}

Measure and fit every region's series, for every grouping key, without the post-fit
stage. This is the pipeline shared by every curve-fit experiment.

The `dataset`/`regions` arguments default to the experiment's own, but can be supplied
explicitly so the GUI can refit live against interactively-positioned regions and noise.
`isfitting=false` (the GUI's Fitting toggle switched off) substitutes [`NoFitting`](@ref)
for the experiment's own model, so the reduced quantities (`x`/`y`) still come through
for the plotted points, but no `curve_fit` call runs and `parameters` comes back empty -
not merely a display toggle, an actual "don't fit" switch.
"""
seriesresults(e::Experiment1D) = seriesresults(e, dataset(e), regions(e))

function seriesresults(e::Experiment1D, ds::Dataset1D, regs; isfitting::Bool=true)
    model = isfitting ? seriesmodel(e) : NoFitting()
    axis = fitaxis(e)
    results = RegionResult[]
    for region in regs
        I = integrate(region, e, ds)
        for (gkey, idx) in groupseries(ds.planes, groupcols(e))
            x = Float64[ds.planes.vars[i][axis] for i in idx]
            y = I[idx]
            perm = sortperm(x)
            x, y, planes = x[perm], y[perm], idx[perm]
            fit = fitseries(model, x, y)
            parameters = OrderedDict{Symbol,Any}(Symbol(n) => p
                                                 for (n, p) in zip(fit.names, fit.params))
            push!(results,
                  RegionResult(region.label, gkey, x, y, planes, parameters, fit.model,
                               fit.converged))
        end
    end
    return results
end

"""
    analyse(e, [dataset, regions]; isfitting=true) -> Vector{RegionResult}

Run the full analysis: reduce, fit, then post-fit. The return type does not depend on the
experiment or on whether anything was fitted, so the GUI's `state[:result]` Observable has
a stable element type even when it starts out empty - which the old
`(; series, summary)` shape did not, `summary` being `nothing` for some experiments and a
`Vector` for others.
"""
analyse(e::Experiment1D) = analyse(e, dataset(e), regions(e))

function analyse(e::Experiment1D, ds::Dataset1D, regs; isfitting::Bool=true)
    results = seriesresults(e, ds, regs; isfitting)
    isfitting || return results
    for label in unique(r.region for r in results)
        postfit!(filter(r -> r.region == label, results), e)
    end
    postfitglobal!(results, e)
    return results
end

# =============================================================================
# shared helpers
# =============================================================================

"""
    Integration(peakppm, noiseppm, ppmwidth)

The integration triple shared with `Exchange1D` (`prob.integration`): a peak position, a
noise position, and a common width, all in ppm (the noise region always has the same
width as the signal region — see [`integrate`](@ref)). Passing one to a top-level
entry point skips the GUI and analyses directly, so a previously-chosen region can be
replayed reproducibly from a script.
"""
const Integration = NamedTuple{(:peakppm, :noiseppm, :ppmwidth)}

regionsfrom(i) = [Region("signal", i.peakppm - i.ppmwidth / 2, i.peakppm + i.ppmwidth / 2)]

"""
    run1d(expt; integration=nothing, call=nothing) -> Vector{RegionResult}

Launch the GUI for `expt` and return the results standing when its window is closed, or -
when an `integration` triple is supplied - skip the GUI and return the analysis for that
region directly.

`call` is the [`AnalysisCall`](@ref) the entry point recorded, carried through to
`summary.txt`. It is unused on the scripted path, which writes no files.

Both paths return the same thing, so what a routine gives back does not depend on how it
was called. [`gui!`](@ref) itself returns the whole GUI state, of which this is one entry;
call it directly where the rest of that state is wanted.
"""
function run1d(expt::Experiment1D; integration=nothing, call=nothing)
    isnothing(integration) && return gui!(expt; call)[:result][]
    d = dataset(expt)
    ds = Dataset1D(d.planes, Float64(integration.noiseppm), d.label, d.sources)
    return analyse(expt, ds, regionsfrom(integration))
end

"""
    defaultregion(dataset; label="signal", width=defaultregionwidth(...)) -> Region

A sensible default single integration region: `width` ppm wide (2% of the spectral
width by default) centred on the tallest peak (by absolute intensity) in the first
plane. Used as the default `regions` for every experiment with a single signal; the GUI
lets the user reposition and resize it, or add further regions.
"""
function defaultregion(dataset::Dataset1D; label="signal",
                       width=defaultregionwidth(first(dataset.planes.traces).δ))
    t = first(dataset.planes.traces)
    peak = t.δ[argmax(abs.(t.y))]
    return Region(label, peak - width / 2, peak + width / 2)
end

# =============================================================================
# implementations - one file per experiment
# =============================================================================

include("expt-relaxation.jl")
include("expt-tract.jl")
include("expt-nutation.jl")
include("expt-diffusion.jl")
include("expt-kinetics.jl")
