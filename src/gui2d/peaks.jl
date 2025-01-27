function Peak(initialposition, label, xradius=0.03, yradius=0.3)
    pos = MaybeVector(initialposition)
    Peak(Observable(pos),
        Observable(label),
        Observable(true), # touched
        Observable(xradius), # xradius
        Observable(yradius), # yradius
        Dict{Symbol, Any}())
end

# generic overlap function - can specialise for different experiments
isoverlapping(peak1, peak2, ::Experiment) = isoverlapping(peak1, peak2)
function isoverlapping(peak1, peak2)
    # Δδs will be a MaybeVector - could be a single element or a list of different shifts
    Δδs = peak1.initialposition[] .- peak2.initialposition[]
    any(Δδ -> begin
        dX = Δδ[1] / (peak1.xradius[] + peak2.xradius[])
        dY = Δδ[2] / (peak1.yradius[] + peak2.yradius[])
        (dX^2 + dY^2) <= 1.0
    end, Δδs)
end