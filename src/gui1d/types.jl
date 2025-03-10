abstract type Experiment end

abstract type FittingModel end
struct NoFitting <: FittingModel end
abstract type ParametricModel <: FittingModel end

struct SpecData
    nmrdata
    dat # data values
    x # x chemical shifts
    y # y axis
    z # z axis (different experiments, may be just [1])
    σ # noise
    xlabel
    ylabel
    zlabels
end

struct Region
    x1
    x2
    label
    integrals
    fitparameters
end

struct Parameter
    label
    value
    uncertainty
end

abstract type VisualisationStrategy end

struct SimpleVisualisation <: VisualisationStrategy end
struct ModelFitVisualisation <: VisualisationStrategy end 
