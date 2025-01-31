function Parameter(label, initialvalue; minvalue=-Inf, maxvalue=Inf, uncertainty=Inf)
    if !isa(initialvalue, MaybeVector)
        initialvalue = MaybeVector(initialvalue)
    end
    value = deepcopy(initialvalue)

    # ensure uncertainty matches size of initialvalue
    if length(uncertainty) != length(value)
        uncertainty = fill(uncertainty, length(value))
    end
    uncertainty = MaybeVector(uncertainty)

    Parameter(label, Observable(value),Observable(uncertainty),
        Observable(initialvalue), Observable(minvalue), Observable(maxvalue))
end

"Return list with value(s) (or if :min or :max passed, their limits)"
function pack!(p, par::Parameter, quantity=:value)
    @debug "Packing parameter $(par.label) ($quantity)" maxlog=10
    x = if quantity == :value
        par.value[]
    elseif quantity == :initial
        par.initialvalue[]
    elseif quantity == :min
        minval = par.minvalue[]
        if length(minval) != length(par.value[])
            fill(minval, length(par.value[]))
        else
            minval
        end
    elseif quantity == :max
        maxval = par.maxvalue[]
        if length(maxval) != length(par.value[])
            fill(maxval, length(par.value[]))
        else
            maxval
        end
    else
        error("Unknown quantity: $quantity")
    end

    append!(p, x)
end

"Unpack a parameter and pop from input vector. Quantity could also be :uncertainty"
function unpack!(v, par::Parameter, quantity=:value)
    n = length(par.value[])
    val = v[1:n]
    deleteat!(v, 1:n)
    # @debug "Unpacking parameter $(par.label) ($quantity): $(first(val))"
        
    if quantity == :value
        par.value[] .= val
    elseif quantity == :uncertainty
        par.uncertainty[] .= val
    end
end

function Base.show(io::IO, p::Parameter)
    print(io, "Parameter($(p.label), value=$(p.value[]), uncertainty=$(p.uncertainty[]))")
end

function Base.show(io::IO, mime::MIME"text/plain", p::Parameter)
    println(io, "Parameter: $(p.label)")
    println(io, "  value: $(p.value[])")
    println(io, "  uncertainty: $(p.uncertainty[])")
    println(io, "  initial value: $(p.initialvalue[])")
    println(io, "  bounds: [$(p.minvalue[]), $(p.maxvalue[])]")
end