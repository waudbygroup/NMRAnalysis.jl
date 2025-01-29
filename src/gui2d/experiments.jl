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
- `specdata`: A SpecData object representing the observed and simulated data plus mask
- `clusters`: An Observable list of clusters of peaks
- `touched`: An Observable list of touched clusters
- `isfitting`: An Observable boolean indicating if real-time fitting is active

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

function setupexptobservables!(expt)
    on(expt.peaks) do _
        mask!(expt)
        cluster!(expt)
        simulate!(expt)
    end
    on(expt.clusters) do _
        checktouched!(expt)
    end
    on(expt.touched, expt.isfitting) do _
        if expt.isfitting[]
            fit!(expt)
        end
    end

    expt
end

function checktouched!(expt)
    @debug "Checking touched clusters"
    touched = Vector{Bool}(length(expt.clusters[]))
    for i in 1:length(expt.clusters[])
        touched[i] = any([expt.peaks[][j].touched[] for j in expt.clusters[i]])
    end
    expt.touched[] = touched
end

function fit!(expt::Experiment)
    @debug "Fitting experiment"
    # iterate over touched clusters and fit
    for i in 1:length(expt.clusters[])
        if expt.touched[][i]
            fit!(expt.clusters[][i], expt)
        end
    end
    notify(expt.peaks)
end

function fit!(cluster::Vector, expt::Experiment)
    @debug "Fitting cluster $cluster"
    # TODO
    # untouch peaks in the cluster
    for i in cluster
        expt.peaks[][i].touched[] = false
    end
end


function simulate!(expt::Experiment)
    @debug "Simulating experiment"
    z = expt.specdata.zfit.val
    for zi in z
        fill!(zi, 0)
    end
    simulate!(z, expt)
    notify(expt.specdata.zfit)
end

function simulate!(z, expt::Experiment)
    @debug "Simulating experiment with z"
    for cluster in expt.clusters[]
        simulate!(z, cluster, expt)
    end
end

function simulate!(z, cluster::Vector, expt::Experiment)
    @debug "Simulating cluster $cluster"
    for i in cluster
        simulate!(z, expt.peaks[][i], expt)
    end
end


function mask!(expt::Experiment)
    @debug "Masking experiment"
    z = expt.specdata.mask.val
    for zi in z
        fill!(zi, false)
    end
    mask!(z, expt)
    notify(expt.specdata.mask)
end

function mask!(z, expt::Experiment)
    @debug "Masking experiment with z"
    for cluster in expt.clusters[]
        mask!(z, cluster, expt)
    end
end

function mask!(z, cluster::Vector, expt::Experiment)
    @debug "Masking experiment with cluster $cluster"
    for i in cluster
        mask!(z, expt.peaks[][i], expt)
    end
end