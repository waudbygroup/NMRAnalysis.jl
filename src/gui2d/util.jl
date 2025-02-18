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

function flatten_with_nan_separator(vectors::Vector{Vector{Point2f}})
    isempty(vectors) && return Point2f[]

    separator = Point2f(NaN, NaN)
    result = reduce(vectors[2:end]; init=vectors[1]) do acc, subvector
        return vcat(acc, [separator], subvector)
    end

    # Remove the trailing separator if it exists
    return length(result) > 0 ? result[1:(end - 1)] : result
end