module NMRAnalysis

using LsqFit
using Measurements
using NMRTools
using Plots
using Reexport
using Statistics

include("viscosity.jl")
include("diffusion.jl")
include("tract.jl")

include("gui2d/GUI2D.jl")
using .GUI2D

using PrecompileTools
include("precompile.jl")

export viscosity
export tract
export diffusion

# @reexport using .GUI: view2d, peakfit2dseries


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
