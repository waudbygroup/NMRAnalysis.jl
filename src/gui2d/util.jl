nearest(A::AbstractArray, t) = findmin(abs.(A .- t))[1]
findnearest(A::AbstractArray, t) = findmin(abs.(A .- t))[2]

function choptitle(title, maxlength=30)
    if length(title) > maxlength
        title[1:maxlength] * "…"
    else
        title
    end
end

function maskellipse!(mask, x, y, x0, y0, xradius, yradius)
    # @debug "masking ellipse at $x0, $y0 with radii $xradius, $yradius" maxlog=10
    x = x .- x0
    y = y' .- y0
    fx = @. yradius^2 * x^2
    fy = @. xradius^2 * y^2
    f = @. fx + fy - xradius^2 * yradius^2
    return mask[f .≤ 0] .= true
end