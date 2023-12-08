module NMRAnalysis

using LsqFit
using Measurements
using NMRTools
using Plots
using Statistics

include("viscosity.jl")
include("diffusion.jl")
include("tract.jl")

using PrecompileTools
include("precompile.jl")

export viscosity
export tract
export diffusion

@info """
NMRAnalysis.jl

Usage:
- set your working directory to a convenient location, e.g.
  cd("/Users/chris/NMR/crick-702/my_experiment_directory")
- call the desired analysis routine

Available analysis routines:
- tract()
- diffusion()

Current working directory: $(pwd())
"""

end
