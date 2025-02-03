module GUI2D

using DelimitedFiles
using GLMakie
using LightGraphs
using LsqFit
using Measurements
using NMRTools
using OrderedCollections


include("util.jl")
include("maybevector.jl")
include("types.jl")
include("parameters.jl")
include("specdata.jl")
include("peaks.jl")
include("experiments.jl")
include("clustering.jl")
include("state.jl")
include("gui.jl")

export MaybeVector

# basic workflow:
# 1. create an Experiment, with associated data
#    - this is handled by the particular experiment type
#    - nmr data are stored as specdata, normalised by ns, rg and concenration, and with a noise value
# 2. the experiment creates a data flow from the peak list (initially empty)
#    to an adjacency matrix, clustering, fitting and simulation of the fitted spectra
#    - this is handled generically for all experiments
#    - a final fitting step can be added on for particular experiment types, e.g. exponential fitting
# 3. create a view state observing the experiment
# 4. add peaks to the experiment to trigger fitting
# 5. modify peaks and fitting will also be triggered

end