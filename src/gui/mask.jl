function getmask(state)
    peakpositions = getpeakpositions(state)
    lift(peakpositions, state[:Xradius], state[:Yradius]) do peakpos, xw, yw
        m = zeros(Bool, size(state[:z][]))
        for pos ∈ peakpos
            maskellipse!(m, state[:x][], state[:y][], pos[1], pos[2], xw, yw)
        end
        m
    end
end

function maskellipse!(mask, δx, δy, δx0, δy0, dx, dy)
    x = δx .- δx0
    y = δy' .- δy0
    f = @. dy^2*x^2 + dx^2*y^2 - dx^2*dy^2
    mask[f .≤ 0] .= true
end
