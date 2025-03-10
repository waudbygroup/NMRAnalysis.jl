module NMRAnalysis

using LsqFit
using Measurements
using NMRTools
using Plots
using Reexport
using Statistics

include("util/util.jl")
using .Util
@reexport using .Util: MaybeVector

include("viscosity.jl")
export viscosity

include("gui1d/GUI1D.jl")
using .GUI1D

# @reexport using .GUI1D: diffusion1d, relaxation1d, modelfit1d

include("gui2d/GUI2D.jl")
using .GUI2D

@reexport using .GUI2D: intensities2d, relaxation2d, recovery2d, modelfit2d # IntensityExperiment
@reexport using .GUI2D: hetnoe2d # HetNOEExperiment
@reexport using .GUI2D: pre2d # PREExperiment

using PrecompileTools
include("precompile.jl")



@info """
NMRAnalysis.jl (v$(pkgversion(NMRAnalysis)))

# 1D Experiment Analysis

- set your working directory to a convenient location, e.g.
  cd("/Users/chris/NMR/crick-702/my_experiment_directory")
- call the desired analysis routine

# 1D experiments

- tract1d()
- diffusion1d()
- relaxation1d()
- modelfit1d()

# 2D experiments

- intensities2d()
- relaxation2d()
- recovery2d()
- modelfit2d()
- hetnoe2d()

Current working directory: $(pwd())
"""

end
