module GUI2D

using CairoMakie
using DelimitedFiles
using GLMakie
using LightGraphs
using LsqFit
using Measurements
using NativeFileDialog
using NMRTools
using OrderedCollections

include("util.jl")
include("maybevector.jl")
include("types.jl")
include("parameters.jl")
include("specdata.jl")
include("peaks.jl")
include("experiments.jl")
include("models.jl")
include("clustering.jl")
include("state.jl")
include("gui.jl")
include("mouse.jl")
include("keyboard.jl")
include("files.jl")
include("visualisation.jl")

export MaybeVector
# export gui!

# export IntensityExperiment
export intensities2d
export relaxation2d
export recovery2d
export modelfit2d

# export HetNOEExperiment
export hetnoe2d

# export CESTExperiment
export cest2d

# export CPMGExperiment
export cpmg2d

# export PREExperiment
export pre2d

end