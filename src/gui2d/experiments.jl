"""
    Experiment

An abstract type representing a type of experimental analysis.

# Implementation Requirements

Concrete subtypes must implement:
- `hasfixedpositions(expt)`: Check if the experiment has fixed peak positions between spectra
- `addpeak!(expt, position)`: Add a peak to the experiment
- `simulate!(z, x, y, peak, expt)`: Simulate a single peak
- `mask!(z, x, y, peak, expt)`: Get the mask for a single peak

Expected fields:
- `peaks`: A list of peaks in the experiment
- `specdata`: A SpecData object representing the simulated data
- `adjacency`: An Observable adjacency matrix for the peaks
- `clusters`: An Observable list of clusters of peaks

# Functions handled by the abstract type

- `nslices(expt)`: Get the number of slices in the experiment
- `npeaks(expt)`: Get the number of peaks in the experiment
- `simulate!(expt)`: Simulate the experiment and update internal specdata
- `fit!(expt)`: Fit the peaks in the experiment

"""

# generic functions
nslices(expt::Experiment) = length(data(expt).z)
npeaks(expt::Experiment) = length(expt.peaks[])

include("expt-relaxation.jl")