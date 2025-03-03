"""
Abstract type representing a 1D NMR experiment.
"""
abstract type Experiment1D end

"""
    SpecData1D

Container for 1D NMR data and metadata.

# Fields
- `nmrdata`: Original NMR data objects
- `x`: Array of x-axis data (chemical shift)
- `y`: Array of y-axis data (intensities)
- `σ`: Noise estimates
- `zlabels`: Labels for slices
- `xplot`: Observable x-axis for display
- `yplot`: Observable y-axis for display
"""
struct SpecData1D
    nmrdata
    x
    y
    σ
    zlabels
    xplot
    yplot
end

"""
    Parameter

Represents a fitting parameter with value, uncertainty, and bounds.

# Fields
- `label`: Name of the parameter
- `value`: Current value (Observable)
- `uncertainty`: Uncertainty in the value (Observable)
- `initialvalue`: Initial value (Observable)
- `minvalue`: Minimum allowed value (Observable)
- `maxvalue`: Maximum allowed value (Observable)
"""
struct Parameter
    label
    value
    uncertainty
    initialvalue
    minvalue
    maxvalue
end

"""
    Region

Represents an integration region in a 1D spectrum.

# Fields
- `label`: Name of the region (Observable)
- `xstart`: Start point on x-axis (Observable)
- `xend`: End point on x-axis (Observable)
- `touched`: Flag indicating if region has been modified (Observable)
- `color`: Color for display (Observable)
- `parameters`: Dictionary of fit parameters
- `postparameters`: Dictionary of derived parameters
- `postfitted`: Flag indicating if post-fit calculations have been done
"""
struct Region
    label::Observable{String}
    xstart::Observable{Float64}
    xend::Observable{Float64}
    touched::Observable{Bool}
    color::Observable{Symbol}
    parameters::OrderedDict{Symbol,Parameter}
    postparameters::OrderedDict{Symbol,Parameter}
    postfitted::Observable{Bool}
end

"""
Abstract type for different fitting models.
"""
abstract type FittingModel end

"""
No fitting, just integration.
"""
struct NoFitting <: FittingModel end

"""
Abstract type for parametric models with fitting functions.
"""
abstract type ParametricModel <: FittingModel end

"""
Abstract type for visualizing data.
"""
abstract type VisualisationStrategy end

"""
Line plot with no offset.
"""
struct LinePlotVisualisation <: VisualisationStrategy end

"""
Stacked plot with vertical offset.
"""
struct StackedPlotVisualisation <: VisualisationStrategy end

"""
Waterfall plot with 3D perspective.
"""
struct WaterfallPlotVisualisation <: VisualisationStrategy end