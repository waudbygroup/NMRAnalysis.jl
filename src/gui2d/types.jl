abstract type Experiment end

struct SpecData
    nmrdata
    x
    y
    z
    σ
    zlabels
    zfit
    mask
end

struct Parameter
    label
    value
    uncertainty
    initialvalue
    minvalue
    maxvalue
end

struct Peak
    parameters::OrderedDict{Symbol, Parameter}
    label::Observable{String}
    touched::Observable{Bool}
    xradius::Observable{Float64}
    yradius::Observable{Float64}
end
