module BatchAnalysis

using CairoMakie
using GLMakie
using NMRTools

include("files.jl")
include("state.jl")
include("linedata.jl")
include("referencing.jl")
include("integration.jl")
include("refgui.jl")
include("integrategui.jl")

export batchanalysis


function batchanalysis(foldername)
    spectra = loadfolder(foldername)
    state = initialisestate(spectra)
    state["foldername"] = foldername
    preparegui!(state)
    ok = referencegui!(state)
    ok || return
    integrategui!(state)

    # return state
end

end