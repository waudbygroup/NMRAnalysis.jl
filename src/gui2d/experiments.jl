"""
    Experiment

An abstract type representing a type of experimental analysis.

# Implementation Requirements

Concrete subtypes must implement:
- `hasfixedpositions(expt)`: Check if the experiment has fixed peak positions between spectra
- `addpeak!(expt, position)`: Add a peak to the experiment
- `simulate!(z, peak, expt)`: Simulate a single peak
- `mask!(z, peak, expt)`: Get the mask for a single peak

Expected fields:
- `peaks`: A list of peaks in the experiment
- `specdata`: A SpecData object representing the simulated data
- `adjacency`: An Observable adjacency matrix for the peaks
- `clusters`: An Observable list of clusters of peaks

# Functions handled by the abstract type

- `nslices(expt)`: Get the number of slices in the experiment
- `npeaks(expt)`: Get the number of peaks in the experiment
- `mask!([z], [peaks], expt)`: Calculate peak masks and update internal specdata
- `simulate!([z], [peaks], expt)`: Simulate the experiment and update internal specdata
- `fit!(expt)`: Fit the peaks in the experiment

"""

# implementations
include("expt-relaxation.jl")

# generic functions
nslices(expt::Experiment) = length(expt.specdata.z)
npeaks(expt::Experiment) = length(expt.peaks[])

function simulate!(expt::Experiment)
    z = expt.specdata.zfit.val
    for zi in z
        fill!(zi, 0)
    end
    simulate!(z, expt)
    notify(expt.specdata.zfit)
end

function simulate!(z, expt::Experiment)
    for cluster in expt.clusters[]
        simulate!(z, cluster, expt)
    end
end

function simulate!(z, cluster, expt::Experiment)
    for i in cluster
        simulate!(z, expt.peaks[][i], expt)
    end
end

function mask!(expt::Experiment)
    @info "Masking experiment"
    z = expt.specdata.mask.val
    for zi in z
        fill!(zi, false)
    end
    mask!(z, expt)
    notify(expt.specdata.mask)
end

function mask!(z, expt::Experiment)
    @info "Masking experiment with z"
    for cluster in expt.clusters[]
        mask!(z, cluster, expt)
    end
end

function mask!(z, cluster::Vector, expt::Experiment)
    @info "Masking experiment with cluster $cluster"
    for i in cluster
        mask!(z, expt.peaks[][i], expt)
    end
end