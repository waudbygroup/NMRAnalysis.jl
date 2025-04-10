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
using ..MCFitting

include("dataset.jl")
include("power.jl")
include("experiments.jl")
include("fitting.jl")
include("Kandkex.jl")
include("state.jl")
include("gui.jl")
include("calibrations.jl")

function r1rho(directory_path="")
    @info "Select a numbered experiment folder or parent directory"
    filenames = select_expts(directory_path; title_filters=["1rho", "1p"])
    isempty(filenames) && return
    return r1rho(filenames)
end

function r1rho(filenames::Vector{String})
    GLMakie.activate!(; focus_on_show=true, title="NMRAnalysis.jl: R1rho fitting")
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