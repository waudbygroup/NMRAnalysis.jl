include("model_amplitudes.jl")
include("model_exp.jl")

updatefitplot(peak, state) =  updatefitplot(peak, state, state[:model])