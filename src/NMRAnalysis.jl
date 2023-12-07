module NMRAnalysis

using LsqFit
using Measurements
using NMRTools
using Plots
using Roots
using Statistics

# include("1d_fitting.jl")
include("tract.jl")

using PrecompileTools
include("precompile.jl")

export analyse_tract

@info """
NMRAnalysis.jl
==============

Usage:
* set your working directory to a convenient location, e.g.
  cd("/Users/chris/NMR/crick-702/my_experiment_directory")
* call the desired analysis routine with the path to the experiment directories

Available analysis routines:
* analyse_tract(trosy_filename, antitrosy_filename)
"""

end
