nearest(A::AbstractArray,t) = A[findmin(abs.(A.-t))[2]]
findnearest(A::AbstractArray,t) = findmin(abs.(A.-t))[2]

function maskellipse!(mask, δH, δN, Hc, Nc, a, b)
    x = δH .- Hc
    y = δN' .- Nc
    f = @. b^2*x^2 + a^2*y^2 - a^2*b^2
    mask[f .≤ 0] .= true
end
