module BatchAnalysis

using NMRTools
using GLMakie

include("files.jl")
include("state.jl")
include("linedata.jl")
include("referencing.jl")
include("refgui.jl")

export batchanalysis


function batchanalysis(foldername)
    spectra = loadfolder(foldername)
    state = initialisestate(spectra)
    preparegui!(state)
    ok = referencegui!(state)

    # return state
end

end