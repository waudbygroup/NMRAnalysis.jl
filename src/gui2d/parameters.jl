function Parameter(label, initialvalue, minvalue=nothing, maxvalue=nothing)
    value = initialvalue
    uncertainty = Inf

    Parameter(label, Observable(value),Observable(uncertainty),
        Observable(initialvalue), Observable(minvalue), Observable(maxvalue))
end

"Return list with value(s) (or if :min or :max passed, their limits)"
function pack(par::Parameter, quantity=:value)
    p = if quantity == :value
        par.value[]
    elseif quantity == :min
        par.minvalue[]
    elseif quantity == :max
        par.maxvalue[]
    elseif quantity == :initial
        par.initialvalue[]
    else
        error("Unknown quantity: $quantity")
    end

    if isarray(p)
        p
    else
        [p]
    end
end

"Unpack a parameter and pop from input vector. Quantity could also be :uncertainty"
function unpack!(v, par::Parameter, quantity=:value)
    val = if isarray(par.value[])
        n = length(par.value[])
        x = v[1:n]
        deleteat!(v, 1:n)
        x
    else
        x = v[1]
        deleteat!(v, 1)
        x
    end
        
    if quantity == :value
        par.value[] = val
    elseif quantity == :uncertainty
        par.uncertainty[] = val
    end
end