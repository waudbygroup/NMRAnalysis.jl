"""
    roiindices(region, trace) -> Vector{Int}

Indices of `trace.δ` falling within `region`. For a zero-width region (or one that
contains no grid points) the single nearest point is returned, giving a peak height.
"""
function roiindices(r::Region, t::Trace)
    idx = findall(δ -> r.lo ≤ δ ≤ r.hi, t.δ)
    isempty(idx) || return idx
    mid = (r.lo + r.hi) / 2
    return [argmin(abs.(t.δ .- mid))]
end

"""
    integrate(region, trace) -> Float64

Summed intensity of `region` in one plane (a height when the region has zero width — the
nearest point is taken). The summed intensity is the convention used by the legacy 1D
routines and by GUI2D's intensity analysis.
"""
integrate(r::Region, t::Trace) = sum(@view t.y[roiindices(r, t)])

"""
    integrate(region, dataset) -> Vector{Measurement{Float64}}

Integrate `region` over every plane of `dataset`, with an uncertainty attached to each.

This is a different noise question from the point-wise spectrum normalisation done in
`tracesfromspec` (dividing by `spec[:noise]`, the whole-spectrum RMS level): the
uncertainty on an *integrated* region depends on its width and on the actual, possibly
non-white noise there, not on a single global scalar. So, following the legacy 1D
routines and `Exchange1D`, it's measured directly: integrate a noise region of the same
width as `region` (centred at `dataset.noisecenter`, which the GUI lets the user
reposition) over every plane, and take the standard deviation of those integrals across
planes. Since the trace intensities are already noise-normalised, this comes out close
to `sqrt(n points)` when the noise is uniform, but tracks the real, locally-measured
noise otherwise.
"""
function integrate(region::Region, ds::Dataset1D)
    w = width(region)
    w == 0 && (w = 0.05)   # a height still needs a nominal window for the noise estimate
    noiseregion = Region("noise", ds.noisecenter - w / 2, ds.noisecenter + w / 2)

    raw = [integrate(region, t) for t in ds.planes.traces]
    noiseintegrals = [integrate(noiseregion, t) for t in ds.planes.traces]
    σ = length(noiseintegrals) > 1 ? std(noiseintegrals) : abs(noiseintegrals[1])
    (σ == 0 || isnan(σ)) && (σ = 1.0)

    return [v ± σ for v in raw]
end
