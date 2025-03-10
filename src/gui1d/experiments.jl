"""
    Experiment

An abstract type representing a 1D experimental analysis.

# Interface Requirements for Concrete Subtypes

Concrete subtypes must implement:

## Required Methods
- `fit(integrals, expt)`: Calculate derived parameters after integration
- `experimentinfo(expt)`: Return formatted text describing experiment
- `regioninfo(expt, idx)`: Return formatted text describing region at index
- `sliceinfo(expt, idx)`: Return description for slice at index

## Required Fields
- `regions::Observable{Vector{Region}}`: List of regions in the experiment
- `specdata::SpecData`: Observed spectral data
- `state::Dict{Symbol,Observable}`: UI and program state

# Generic Functionality (Provided by Abstract Base)

The following methods are implemented for all experiment types:
    
## Measurement and Access
- `nslices(expt)`: Number of spectra in the experiment
- `nregions(expt)`: Number of regions currently in the experiment
- `visualisationtype(expt)`: Returns the visualization strategy for the experiment (can be overridden)
    
## Region Manipulation
- `addregion!(expt, x1, x2, [label::String])`: Add an integration region at specified position
- `deleteregion!(expt, idx)`: Delete region at specified index
- `deleteallregions!(expt)`: Remove all regions from experiment
"""

nslices(expt::Experiment) = length(expt.specdata.z)
nregions(expt::Experiment) = length(expt.regions[])
visualisationtype(::Experiment) = SimpleVisualisation()
setupexptobservables!(::Experiment) = nothing
fit(_, ::Experiment) = Dict{Symbol, Parameter}()

function integrate(x1, x2, expt::Experiment)
    xl = min(x1, x2)
    xh = max(x1, x2)
    sd = expt.specdata
    map(1:nslices(expt)) do i
        sum(sd.y[i][findall(x -> xl <= x <= xh, sd.x[i]),:], dims=1)
    end
end

function addregion!(expt, x1, x2, label="")
    push!(expt.regions[], Region(x1, x2, label, expt))
    expt.state[:current_peak_idx][] = length(expt.regions[])
end

function deleteregion!(expt, idx)
    # update current peak index if necessary
    if expt.state[:current_peak_idx][] == idx
        expt.state[:current_peak_idx][] = 0
    elseif expt.state[:current_peak_idx][] > idx
        expt.state[:current_peak_idx][] -= 1
    end
    # delete row idx from regions matrix
    expt.regions[] = vcat(expt.regions[][1:idx-1, :], expt.regions[][idx+1:end, :])
end

function deleteallregions!(expt)
    expt.state[:current_peak_idx][] = 0
    expt.regions[] = zeros(0, 2)
end