using NMRAnalysis
using Test
using SafeTestsets

@safetestset "diffusion test" begin
    include("diffusion_test.jl")
end
