"Find non-overlapping clusters for fitting, return list of clusters (e.g. [[1], [2, 3]])"
function makeclusters(state)
    peakradiusX = state[:Xradius][]
    peakradiusY = state[:Yradius][]
    
    @debug "Forming non-overlapping clusters of peaks for fitting" peakradiusX peakradiusY

    # function to test for overlap between peaks i and j
    # overlaps(i, j) = (abs(peaks[i,1]-peaks[j,1]) < peakradiusH) && (abs(peaks[i,2]-peaks[j,2]) < peakradiusN)
    function overlaps(i, j)
        Δδs = state[:peaks][][i].initial_position .- state[:peaks][][j].initial_position
        any(map(Δδs) do Δδ
            dX = Δδ[1] / (2*peakradiusX)
            dY = Δδ[2] / (2*peakradiusY)
            (dX^2 + dY^2) <= 1.0
        end)
        
    end

    npeaks = length(state[:peaks][])
    adjacencymatrix = [overlaps(i,j) for i=1:npeaks, j=1:npeaks]
    clusters = connected_components(SimpleGraph(adjacencymatrix))
    
    @debug "$(length(clusters)) clusters of peaks found" clusters
    return clusters
end