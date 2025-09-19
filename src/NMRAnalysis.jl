module NMRAnalysis

using LsqFit
using Measurements
using NativeFileDialog
using NMRTools
using Plots
using REPL.TerminalMenus
using Reexport
using Statistics

include("fileselection.jl")
include("viscosity.jl")
include("diffusion.jl")
include("relaxation.jl")
include("tract.jl")

include("gui2d/GUI2D.jl")
using .GUI2D

using PrecompileTools
include("precompile.jl")

export viscosity
export diffusion
export relaxation
export tract

include("R1rho/R1rho.jl")
using .R1rho

@reexport using .GUI2D: MaybeVector
@reexport using .GUI2D: intensities2d, relaxation2d, recovery2d, modelfit2d # IntensityExperiment
@reexport using .GUI2D: hetnoe2d # HetNOEExperiment
@reexport using .GUI2D: cest2d # CESTExperiment
@reexport using .GUI2D: cpmg2d # CPMGExperiment
@reexport using .GUI2D: pre2d # PREExperiment

@reexport using .R1rho: r1rho, setupR1rhopowers

@info """
NMRAnalysis.jl (v$(pkgversion(NMRAnalysis)))

# 1D Experiment Analysis

- set your working directory to a convenient location, e.g.
  cd("/Users/chris/NMR/crick-702/my_experiment_directory")
- call the desired analysis routine

Available analysis routines:
- tract()
- diffusion()
- r1rho()

# 2D Experiment Analysis

Available interfaces:
- intensities2d()
- relaxation2d()
- recovery2d()
- modelfit2d()
- hetnoe2d()
- cest2d()
- cpmg2d()

Current working directory: $(pwd())
"""

end
