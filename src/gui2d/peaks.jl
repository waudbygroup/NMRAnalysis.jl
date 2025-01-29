function Peak(initialposition, label, xradius=0.03, yradius=0.3)
    pos = MaybeVector(initialposition)
    Peak(Observable(pos),
        Observable(label),
        Observable(true), # touched
        Observable(xradius), # xradius
        Observable(yradius), # yradius
        OrderedDict{Symbol, Any}())
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

"Create vector of parameters (or if :min or :max passed, their limits)"
function pack(peak::Peak, quantity=:value)
    # 1. positions
    # TODO - need to rework this to store fitted positions
    xy = peak.initialposition[]
    x = [pos[1] for pos in xy]
    y = [pos[2] for pos in xy]
    if quantity == :min
        x .-= peak.xradius[]
        y .-= peak.yradius[]
    elseif quantity == :max
        x .+= peak.xradius[]
        y .+= peak.yradius[]
    end
    p = [x; y]

    # 2. parameter dictionary
    d = peak.parameters
    for par in values(d)
        append!(p, pack(par, quantity))
    end

    p
end

"Unpack a parameter and pop from input vector. Quantity could also be :uncertainty"
function unpack!(v, peak::Peak, quantity=:value)
    # 1. positions
    # TODO - need to rework this to store fitted positions
    xy = peak.initialposition[]
    n = length(x)
    x = v[1:n]
    y = v[n+1:2n]
    deleteat!(v, 1:2n)

    for i=1:n
        xy[i] = Point2f(x[i], y[i])
    end
    peak.initialposition[] = [x y]

    # 2. parameter dictionary
    d = peak.parameters
    for par in values(d)
        unpack!(v, par, quantity)
    end
end