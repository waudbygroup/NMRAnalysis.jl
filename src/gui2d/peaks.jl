function Peak(initialposition, label, xradius=0.03, yradius=0.3)
    Peak(Observable(initialposition),
        Observable(label),
        Observable(true), # touched
        Observable(xradius), # xradius
        Observable(yradius), # yradius
        DefaultDict{Symbol, Any}(nothing))
end

# generic overlap function - can specialise for different experiments
isoverlapping(peak1, peak2, expt::Experiment) = isoverlapping(peak1, peak2, typeof(expt))
isoverlapping(peak1, peak2) = isoverlapping(peak1, peak2, ::Experiment)
function isoverlapping(peak1, peak2, ::Experiment)
    # Δδs will be a MaybeVector - could be a single element or a list of different shifts
    Δδs = peak1.initialposition[] .- peak2.initialposition[]
    any(Δδ -> begin
        dX = Δδ[1] / (peak1.xradius[] + peak2.xradius[])
        dY = Δδ[2] / (peak1.yradius[] + peak2.yradius[])
        (dX^2 + dY^2) <= 1.0
    end, Δδs)
end