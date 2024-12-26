"Find non-overlapping clusters for fitting, return list of clusters (e.g. [[1], [2, 3]])"
function makeclusters(peaks) #, setup)
    peakradiusX = 0.05 #setup.options.peakradiusH
    peakradiusY = 0.5 #setup.options.peakradiusN
    
    @debug "Forming non-overlapping clusters of peaks for fitting" peakradiusH peakradiusN

    # function to test for overlap between peaks i and j
    # overlaps(i, j) = (abs(peaks[i,1]-peaks[j,1]) < peakradiusH) && (abs(peaks[i,2]-peaks[j,2]) < peakradiusN)
    function overlaps(i, j)
        Δδ = peaks[][i].position - peaks[][j].position
        dX = Δδ[1] / (2*peakradiusX)
        dY = Δδ[2] / (2*peakradiusY)
        (dX^2 + dY^2) <= 1.0
    end

    npeaks = length(peaks[])
    adjacencymatrix = [overlaps(i,j) for i=1:npeaks, j=1:npeaks]
    clusters = connected_components(SimpleGraph(adjacencymatrix))
    # e.g. clusters = [[1], [2, 3]]
    
    @debug "$(length(clusters)) clusters of peaks found" clusters
    return clusters
end