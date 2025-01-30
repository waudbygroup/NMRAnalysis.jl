function Parameter(label, initialvalue, minvalue=nothing, maxvalue=nothing)
    value = initialvalue
    uncertainty = Inf

    Parameter(label, Observable(value),Observable(uncertainty),
        Observable(initialvalue), Observable(minvalue), Observable(maxvalue))
end

"Return list with value(s) (or if :min or :max passed, their limits)"
function pack!(p, par::Parameter, quantity=:value)
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
        
    if quantity == :value
        par.value[] .= val
    elseif quantity == :uncertainty
        par.uncertainty[] .= val
    end
end