module R1rho

using CairoMakie
using GLMakie
using LinearAlgebra
using LsqFit
using Measurements
using NMRTools
using Statistics

export r1rho

include("dataset.jl")
include("power.jl")
include("experiments.jl")
include("fitting.jl")
include("state.jl")
include("gui.jl")

function r1rho(filenames)
    GLMakie.activate!()
    dataset = processexperiments(filenames)
    state = initialisestate(dataset)
    state[:filenames] = filenames
    gui!(state)
end

end