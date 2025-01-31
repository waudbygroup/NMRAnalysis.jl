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
        @debug "Peaks changed"
        mask!(expt)
        cluster!(expt)
        simulate!(expt)
    end
    on(expt.clusters) do _
        @debug "Clusters changed"
        checktouched!(expt)
    end
    on(expt.touched) do _
        @debug "Cluster touched"
        # touch peaks inside touched clusters
        for i in 1:length(expt.touched[])
            if expt.touched[][i]
                for j in expt.clusters[][i]
                    expt.peaks[][j].touched[] = true
                end
            end
        end
        if expt.isfitting[]
            fit!(expt)
        end
    end
    on(expt.isfitting) do _
        @debug "Fitting changed"
        if expt.isfitting[]
            fit!(expt)
        end
    end

    expt
end

function deletepeak!(expt, idx)
    # touch other peaks in the same cluster before deleting and notifying
    clusteridx = findfirst(i -> idx in i, expt.clusters[])
    cluster = expt.clusters[][clusteridx]
    for i in cluster
        expt.peaks[][i].touched[] = true
    end
    deleteat!(expt.peaks[], idx)
    notify(expt.peaks)
end

function checktouched!(expt)
    @debug "Checking touched clusters" maxlog=10
    touched = map(expt.clusters[]) do cluster
        any([expt.peaks[][j].touched[] for j in cluster])
    end
    expt.touched[] = touched
end

function fit!(expt::Experiment)
    @debug "Fitting experiment" maxlog=10
    anythingchanged = false
    # iterate over touched clusters and fit
    for i in 1:length(expt.clusters[])
        if expt.peaks[][i].touched[]
            fit!(expt.clusters[][i], expt)
            postfit!(expt.clusters[][i], expt)
            anythingchanged = true
        end
    end
    if anythingchanged
        @debug "Fit finished - notifying peaks" maxlog=10
        postfitglobal!(expt)
        notify(expt.peaks)
    end
end

# placeholder functions for additional fitting following spectrum fit
function postfit!(cluster::Vector{Int}, expt::Experiment)
    @debug "Post-fitting cluster $cluster" maxlog=10
    for i in cluster
        postfit!(expt.peaks[][i], expt)
    end
end

"Additional fitting of peak following spectrum fit - defaults to no action"
postfit!(peak::Peak, expt::Experiment) = nothing

"Global fitting of entire experiment following spectrum fit - defaults to no action"
postfitglobal!(expt::Experiment) = nothing

function fit!(cluster::Vector{Int}, expt::Experiment)
    @debug "Fitting cluster $cluster" #maxlog=10
    peaks = [expt.peaks[][i] for i in cluster]

    # TODO - adjust to work with smaller area of spectra

    # initial parameters
    p0 = pack(peaks, :initial)
    pmin = pack(peaks, :min)
    pmax = pack(peaks, :max)

    m = mask(cluster, expt)
    zobs = reduce(vcat, vec.(expt.specdata.z))[m]
    
    zsim = [similar(zi) for zi in expt.specdata.z]
    zsimm = similar(zobs)
    # create residual function
    function resid(p)
        @debug "resid (start)" maxlog=10
        unpack!(copy(p), peaks, :value)
        for i=1:nslices(expt)
            fill!(zsim[i], 0.0)
        end
        simulate!(zsim, cluster, expt)
        zsimm .= reduce(vcat, vec.(zsim))[m]
        zobs - zsimm
    end
    @debug "running fit" maxlog=10
    sol = LsqFit.lmfit(resid, p0, Float64[])
    pfit = coef(sol)
    perr = stderror(sol)
    @debug "fit complete" pfit maxlog=10
    unpack!(pfit, peaks, :value)
    unpack!(perr, peaks, :uncertainty)

    # untouch peaks in the cluster
    for peak in peaks
        peak.touched.val = false
    end
    # N.B. the update to peaks[] will be notified in parent function once all fitting is complete
end


function simulate!(expt::Experiment)
    @debug "Simulating experiment" maxlog=10
    z = expt.specdata.zfit.val
    for zi in z
        fill!(zi, 0)
    end
    simulate!(z, expt)
    notify(expt.specdata.zfit)
end

function simulate!(z, expt::Experiment)
    @debug "Simulating experiment with z" maxlog=10
    for cluster in expt.clusters[]
        simulate!(z, cluster, expt)
    end
end

function simulate!(z, cluster::Vector{Int}, expt::Experiment)
    @debug "Simulating cluster $cluster" maxlog=10
    for i in cluster
        simulate!(z, expt.peaks[][i], expt)
    end
end


function mask!(expt::Experiment)
    @debug "Masking experiment" maxlog=10
    z = expt.specdata.mask.val
    for zi in z
        fill!(zi, false)
    end
    mask!(z, expt)
    notify(expt.specdata.mask)
end

function mask!(z, expt::Experiment)
    @debug "Masking experiment with z" maxlog=10
    for peak in expt.peaks[]
        mask!(z, peak, expt)
    end
end

function mask(cluster::Vector{Int}, expt::Experiment)
    @debug "Generating cluster mask" maxlog=10
    z = [similar(zi) for zi in expt.specdata.mask[]]
    for zi in z
        fill!(zi, false)
    end
    for i in cluster
        mask!(z, expt.peaks[][i], expt)
    end
    reduce(vcat, vec.(z))
end

function Base.show(io::IO, expt::Experiment)
    print(io, "$(typeof(expt))($(npeaks(expt)) peaks, $(nslices(expt)) slices)")
end

function Base.show(io::IO, mime::MIME"text/plain", expt::Experiment)
    println(io, "$(typeof(expt))")
    println(io, "  $(npeaks(expt)) peaks")
    println(io, "  $(nslices(expt)) slices")
    println(io, "  $(length(expt.clusters[])) clusters")
    println(io, "  fitting: $(expt.isfitting[])")
end