module R1rho

using CairoMakie
using Distributions: cdf, FDist
using GLMakie
using LinearAlgebra
using LsqFit
using MonteCarloMeasurements: ±, pmean, pstd, register_primitive
using NMRTools
using Printf
using Random
using Statistics

export r1rho, setupR1rhopowers
using ..NMRAnalysis: select_expts

include("dataset.jl")
include("power.jl")
include("experiments.jl")
include("fitting.jl")
include("Kandkex.jl")
include("state.jl")
include("gui.jl")
include("calibrations.jl")

"""
    r1rho(directory_path=""; scalefactor=:automatic)
    r1rho(filenames::Vector{String}; scalefactor=:automatic)

Launch the R1rho analysis GUI for a given directory or experiment folder.

- `directory_path`: Path to the experiment folder or parent directory. If not specified, a dialog will prompt for selection.
- `filenames`: Vector of experiment files to analyze.
- `scalefactor`: Adjusts the size of the display window. By default, this is `2` for high-resolution displays and `1` for low-resolution displays. Omit or set to `:automatic` to use the default, or provide a numeric value to override.

# Example

```julia
# analyse experiments 104 and 105
r1rho(["/Users/chris/git/R1rho.jl/data/104", "/Users/chris/git/R1rho.jl/data/105"])

# open a dialog to select a directory, and scale the display size by 1.5
r1rho(scalefactor=1.5)
```
"""
function r1rho(directory_path=""; scalefactor=:automatic)
    @info "Select a numbered experiment folder or parent directory"
    filenames = select_expts(directory_path; title_filters=["1rho", "1p"])
    isempty(filenames) && return
    return r1rho(filenames; scalefactor=scalefactor)
end

function r1rho(filenames::Vector{String}; scalefactor=:automatic)
    if scalefactor == :automatic
        GLMakie.activate!(; focus_on_show=true, title="NMRAnalysis.jl: R1rho fitting")
    elseif scalefactor isa Number
        GLMakie.activate!(; focus_on_show=true, title="NMRAnalysis.jl: R1rho fitting",
                          scalefactor=scalefactor)
    else
        @error "scalefactor must be :auto or a number"
    end
    dataset = processexperiments(filenames)
    state = initialisestate(dataset)
    state[:filenames] = filenames
    return gui!(state)
end

function __init__()
    # Register primitives here - this gets called when the module is loaded
    # register the function with MonteCarloMeasurements as a primative
    return register_primitive(safesqrt)
end

end