module NMRAnalysis

using LsqFit
using Measurements
using NativeFileDialog
using NMRTools
using Plots
using REPL.TerminalMenus
using Statistics

include("fileselection.jl")
include("viscosity.jl")
include("diffusion.jl")
include("tract.jl")

using PrecompileTools
include("precompile.jl")

export viscosity
export tract
export diffusion

include("R1rho/R1rho.jl")
using .R1rho
export r1rho


@info """
NMRAnalysis.jl

Usage:
- set your working directory to a convenient location, e.g.
  cd("/Users/chris/NMR/crick-702/my_experiment_directory")
- call the desired analysis routine

Available analysis routines:
- tract()
- diffusion()
- r1rho()

Current working directory: $(pwd())
"""

end
