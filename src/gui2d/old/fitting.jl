function fit!(state)
    # are there any touched peaks to be fitted?
    anytouched = any([peak.touched for peak ∈ state[:peaks][]])

    # iterate over untouched peaks
    # we will fit clusters of peaks together
    # we need to populate the dictionary sim_parameters within each peak,
    # include the position, amplitude, and linewidths
    # these can be used by the simZ function to generate the simulated spectrum
    if anytouched
        for cluster in state[:clusters]
            # skip over untouched clusters
            any([state[:peaks][][peakindex].touched for peakindex ∈ cluster]) || continue
            fitcluster!(state, cluster, state[:model])
            for peakindex in cluster
                state[:peaks][][peakindex].touched = false
            end
        end
        @debug "fit! done fitting"
    end

    simZ!(state)
    state[:zfit][] = state[:specdata][:zfit][state[:slice][]]
    state[:fitpositions][] = getfitpositions(state)
    @debug "fit! done notifying zfit"
end
