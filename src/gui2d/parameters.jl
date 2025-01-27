function Parameter(label, initialvalue, minvalue=nothing, maxvalue=nothing)
    value = initialvalue
    uncertainty = Inf

    Parameter(label, Observable(value),Observable(uncertainty),
        Observable(initialvalue), Observable(minvalue), Observable(maxvalue))
end

# TODO pack/unpack functions