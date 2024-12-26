module GUI

using DelimitedFiles
using GLMakie
using LightGraphs
using NMRTools

export view2d
export peakfit2dseries

include("util.jl")
include("types.jl")
include("peaks.jl")
include("clustering.jl")
include("specdata.jl")
include("view2d.jl")
include("peakfit2dseries.jl")
include("simpeak.jl")
# include("fit2d.jl")
# include("mouse.jl")
# include("keyboard.jl")
# include("mask.jl")
# include("models.jl")
# include("sim.jl")

end