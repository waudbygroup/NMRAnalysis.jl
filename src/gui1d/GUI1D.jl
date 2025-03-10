module GUI1D

using CairoMakie
using DelimitedFiles
using GLMakie
using LsqFit
using Measurements
using NativeFileDialog
using NMRTools
using OrderedCollections
using Statistics

using ..Util

# include("util.jl")
include("types.jl")
include("region.jl")
include("models.jl")
include("experiments.jl")
include("state.jl")


# Include experiment-specific implementations

include("expt-curvefit.jl")
# Export functions that will be available to the user
# export diffusion1d
# export tract1d
# export relaxation1d
# export modelfit1d

end