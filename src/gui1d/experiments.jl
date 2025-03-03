"""
    Experiment1D

An abstract type representing a type of 1D experimental analysis.

# Implementation Requirements

Concrete subtypes must implement:
- `integrate!(region::Region, expt)`: Calculate integral for a region
- `fit!(region::Region, expt)`: Fit the region data
- `postfit!(region::Region, expt)`: Perform post-fitting calculations

Expected fields:
- `specdata`: A SpecData1D object representing the observed data
- `regions`: An Observable list of integration regions
- `isfitting`: An Observable boolean indicating if real-time fitting is active
- `state`: An Observable dictionary of state variables
- `colors`: A vector of colors for regions
"""

"""
    nslices(expt::Experiment1D)

Number of spectra in the experiment.
"""
nslices(expt::Experiment1D) = length(expt.specdata.y)

"""
    nregions(expt::Experiment1D) 

Number of regions currently in the experiment.
"""
nregions(expt::Experiment1D) = length(expt.regions[])

"""
    setupexptobservables!(expt)

Set up observables and callbacks for the experiment.
"""
function setupexptobservables!(expt)
    on(expt.regions) do _
        @debug "Regions changed"
        checktouched!(expt)
    end
    
    on(expt.isfitting) do _
        @debug "Fitting status changed"
        if expt.isfitting[]
            fit!(expt)
        end
    end
    
    return expt
end

"""
    checktouched!(expt)

Update which regions have been modified.
"""
function checktouched!(expt::Experiment1D)
    @debug "Checking touched regions" 
    
    # Check if any region is touched and needs fitting
    if expt.isfitting[]
        anytouched = any(region -> region.touched[], expt.regions[])
        if anytouched
            fit!(expt)
        end
    end
end

"""
    preparespecdata(inputfilenames, ::Type{E}) where E <: Experiment1D

Load NMR data and prepare the SpecData1D object.
"""
function preparespecdata(inputfilenames, ::Type{E}) where E <: Experiment1D
    @debug "Preparing spec data for experiment: $inputfilenames"
    
    # Handle single file case
    if inputfilenames isa String
        return loadsingleexperiment(inputfilenames, E)
    end
    
    # Handle multiple files case
    if inputfilenames isa Vector{String}
        return loadmultipleexperiments(inputfilenames, E)
    end
    
    error("Unsupported input format")
end

"""
    loadsingleexperiment(inputfilename, ::Type{E}) where E <: Experiment1D

Load a single NMR experiment and prepare SpecData1D.
"""
function loadsingleexperiment(inputfilename, ::Type{E}) where E <: Experiment1D
    @debug "Loading single experiment: $inputfilename"
    
    # Load NMR data
    nmrdata = loadnmr(inputfilename)
    
    # Extract data
    if ndims(nmrdata) == 2  # Pseudo-2D
        # Get axis data
        x = data(nmrdata, F1Dim)
        
        # Get intensity data for each slice
        y = [data(nmrdata, i) for i in 1:size(nmrdata, 2)]
        
        # Estimate noise for each slice
        σ = [nmrdata[:noise] for _ in 1:length(y)]
        
        # Create slice labels
        zlabels = ["Slice $i" for i in 1:length(y)]
    else  # 1D
        # Get axis data
        x = [data(nmrdata, F1Dim)]
        
        # Get intensity data
        y = [data(nmrdata)]
        
        # Estimate noise
        σ = [nmrdata[:noise]]
        
        # Create slice label
        zlabels = ["Spectrum"]
    end
    
    # Create observables for plotting
    xplot = Observable(x[1])
    yplot = Observable(y[1])
    
    return SpecData1D(nmrdata, x, y, σ, zlabels, xplot, yplot)
end

"""
    loadmultipleexperiments(inputfilenames, ::Type{E}) where E <: Experiment1D

Load multiple NMR experiments and prepare SpecData1D.
"""
function loadmultipleexperiments(inputfilenames, ::Type{E}) where E <: Experiment1D
    @debug "Loading multiple experiments: $inputfilenames"
    
    # Load all NMR data
    nmrdata_list = map(loadnmr, inputfilenames)
    
    # Extract data
    x = map(nmr -> data(nmr, F1Dim), nmrdata_list)
    y = map(data, nmrdata_list)
    σ = map(nmr -> nmr[:noise], nmrdata_list)
    
    # Create slice labels from filenames
    zlabels = map(nmr -> basename(nmr[:filename]), nmrdata_list)
    
    # Create observables for plotting
    xplot = Observable(x[1])
    yplot = Observable(y[1])
    
    return SpecData1D(nmrdata_list, x, y, σ, zlabels, xplot, yplot)
end

"""
    Base.show(io::IO, expt::Experiment1D)

Display basic information about the experiment.
"""
function Base.show(io::IO, expt::Experiment1D)
    print(io, "$(typeof(expt))($(nregions(expt)) regions, $(nslices(expt)) slices)")
end

"""
    Base.show(io::IO, mime::MIME"text/plain", expt::Experiment1D)

Display detailed information about the experiment.
"""
function Base.show(io::IO, mime::MIME"text/plain", expt::Experiment1D)
    println(io, "$(typeof(expt))")
    println(io, "  $(nregions(expt)) regions")
    println(io, "  $(nslices(expt)) slices")
    println(io, "  fitting: $(expt.isfitting[])")
end

"""
    renameregion!(expt, state, initiator)

Start the process of renaming a region.
"""
function renameregion!(expt, state, initiator)
    @debug "Renaming region"
    
    if initiator == :keyboard
        state[:mode][] = :renamingstart
    else
        state[:mode][] = :renaming
    end
    
    state[:oldlabel][] = state[:current_region][].label[]
    state[:current_region][].label[] = "‸"
    notify(expt.regions)
end

"""
    experimentinfo(expt::Experiment1D)

Return formatted text describing the experiment.
"""
function experimentinfo(expt::Experiment1D)
    return "Analysis type: $(typeof(expt))\n" *
           "Number of regions: $(nregions(expt))\n" *
           "Number of slices: $(nslices(expt))"
end

"""
    slicelabel(expt::Experiment1D, idx)

Return descriptive text for slice idx.
"""
function slicelabel(expt::Experiment1D, idx)
    if idx > 0 && idx <= length(expt.specdata.zlabels)
        return "$(expt.specdata.zlabels[idx]) ($idx of $(nslices(expt)))"
    else
        return "Invalid slice index: $idx"
    end
end

"""
    regioninfotext(expt::Experiment1D, idx)

Return formatted text describing region idx.
"""
function regioninfotext(expt::Experiment1D, idx)
    if idx == 0
        return "No region selected"
    end
    
    region = expt.regions[][idx]
    
    if !region.postfitted[]
        return "Region: $(region.label[])\nNot fitted"
    end
    
    # Common region information
    info = [
        "Region: $(region.label[])",
        "Range: $(round(region.xstart[], digits=3)) - $(round(region.xend[], digits=3)) ppm",
        "Width: $(round(width(region), digits=3)) ppm"
    ]
    
    # Add model-specific parameters if any
    if !isempty(region.postparameters)
        push!(info, "")
        push!(info, "Fitted parameters:")
        for (_, param) in region.postparameters
            push!(info, "$(param.label): $(param.value[][1]) ± $(param.uncertainty[][1])")
        end
    end
    
    return join(info, "\n")
end

"""
    update_specdata_display!(expt::Experiment1D)

Update the displayed spectrum data.
"""
function update_specdata_display!(expt::Experiment1D)
    # Get current slice
    slice_idx = expt.state[][:current_slice][]
    
    # Update x and y data
    if slice_idx > 0 && slice_idx <= nslices(expt)
        expt.specdata.xplot[] = expt.specdata.x[slice_idx]
        expt.specdata.yplot[] = expt.specdata.y[slice_idx]
    end
end