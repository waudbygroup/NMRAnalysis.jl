module NMRAnalysis

using LsqFit
using Measurements
using NMRTools
using Plots
using Roots
using Statistics

include("diffusion.jl")
include("tract.jl")

using PrecompileTools
include("precompile.jl")

export analyse_tract
export analyse_diffusion

@info """
NMRAnalysis.jl

Usage:
- set your working directory to a convenient location, e.g.
  cd("/Users/chris/NMR/crick-702/my_experiment_directory")
- call the desired analysis routine

Available analysis routines:
- analyse_tract()
- analyse_diffusion()

Current working directory: $(pwd())
"""

end
