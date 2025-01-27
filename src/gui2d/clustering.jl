function makeadjacency(peaks, expt)
    npeaks = length(peaks)
    adjacency = zeros(npeaks, npeaks)

    for i in 1:npeaks
        for j in i+1:npeaks
            adjacency[i, j] = isoverlapping(peaks[i], peaks[j], expt)
            adjacency[j, i] = adjacency[i, j]
        end
    end

    adjacency
end

"Find non-overlapping clusters for fitting, return list of clusters (e.g. [[1], [2, 3]])"
function findclusters(adjacencymatrix)
    clusters = connected_components(SimpleGraph(adjacencymatrix))
    @debug "$(length(clusters)) clusters of peaks found" clusters

    clusters
end
