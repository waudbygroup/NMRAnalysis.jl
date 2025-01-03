module R1rho

using CairoMakie
using GLMakie
using LinearAlgebra
using LsqFit
using Measurements
using NMRTools
using Statistics

export r1rho
using ..NMRAnalysis: select_expts

include("dataset.jl")
include("power.jl")
include("experiments.jl")
include("fitting.jl")
include("state.jl")
include("gui.jl")

function r1rho(directory_path="")
    @info "Select a numbered experiment folder or parent directory"
    filenames = select_expts(directory_path; title_filters=["1rho","1p"])
    isempty(filenames) && return
    r1rho(filenames)
end

function r1rho(filenames::Vector{String})
    GLMakie.activate!(;focus_on_show=true, title="NMRAnalysis.jl: R1rho fitting")
    dataset = processexperiments(filenames)
    state = initialisestate(dataset)
    state[:filenames] = filenames
    gui!(state)
end

end