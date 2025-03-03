"""
    Region(xstart, xend, label="", color=:blue)

Create a new integration region with the given start and end points.
"""
function Region(xstart, xend, label="", color=:blue)
    # Ensure start is less than end
    xmin, xmax = minmax(xstart, xend)
    
    # Create observable values
    label_obs = Observable(label)
    xstart_obs = Observable(xmin)
    xend_obs = Observable(xmax)
    touched_obs = Observable(true)
    color_obs = Observable(color)
    
    # Create parameter dictionaries
    params = OrderedDict{Symbol,Parameter}()
    postparams = OrderedDict{Symbol,Parameter}()
    
    # Default parameter for integral
    params[:integral] = Parameter("Integral", 0.0)
    
    return Region(
        label_obs,
        xstart_obs,
        xend_obs,
        touched_obs,
        color_obs,
        params,
        postparams,
        Observable(false)
    )
end

"""
    addregion!(expt::Experiment1D, position)

Add a new integration region to the experiment, centered at `position`
with a default width.
"""
function addregion!(expt::Experiment1D, position::Float64)
    # Default width is 0.2 ppm
    halfwidth = 0.1
    xstart = position - halfwidth
    xend = position + halfwidth
    
    expt.state[][:total_regions][] += 1
    label = "R$(expt.state[][:total_regions][])"
    
    # Get next color from cycle
    color_idx = (expt.state[][:total_regions][] - 1) % length(expt.colors) + 1
    color = expt.colors[color_idx]
    
    # Create the region
    newregion = Region(xstart, xend, label, color)
    
    # Add parameters based on experiment type
    setup_region_parameters!(newregion, expt)
    
    # Add to the regions list
    push!(expt.regions[], newregion)
    notify(expt.regions)
    
    # Set as current region
    expt.state[][:current_region_idx][] = length(expt.regions[])
    
    return newregion
end

"""
    moveregion!(expt, idx, delta)

Move region `idx` by `delta` ppm. Updates fit automatically.
"""
function moveregion!(expt, idx, delta)
    region = expt.regions[][idx]
    region.touched[] = true
    region.xstart[] += delta
    region.xend[] += delta
    notify(expt.regions)
end

"""
    resizeregion!(expt, idx, isstart, delta)

Resize region `idx` by moving either the start or end point by `delta`.
"""
function resizeregion!(expt, idx, isstart, delta)
    region = expt.regions[][idx]
    region.touched[] = true
    
    if isstart
        region.xstart[] += delta
    else
        region.xend[] += delta
    end
    
    # Ensure start is always less than end
    if region.xstart[] > region.xend[]
        xstart = region.xend[]
        xend = region.xstart[]
        region.xstart[] = xstart
        region.xend[] = xend
    end
    
    notify(expt.regions)
end

"""
    deleteregion!(expt, idx)

Delete region `idx`. Updates fit automatically.
"""
function deleteregion!(expt, idx)
    region = expt.regions[][idx]
    region.touched[] = true
    deleteat!(expt.regions[], idx)
    
    # Update current region index
    if expt.state[][:current_region_idx][] >= idx
        expt.state[][:current_region_idx][] = max(0, expt.state[][:current_region_idx][] - 1)
    end
    
    notify(expt.regions)
end

"""
    deleteallregions!(expt)

Remove all regions from the experiment.
"""
function deleteallregions!(expt)
    expt.state[][:current_region_idx][] = 0
    # Delete any existing regions
    for i in length(expt.regions[]):-1:1
        deleteregion!(expt, i)
    end
    expt.state[][:total_regions][] = 0
end

"""
    integrate!(expt::Experiment1D)

Calculate integrals for all regions in the experiment.
"""
function integrate!(expt::Experiment1D)
    for region in expt.regions[]
        integrate!(region, expt)
    end
    
    notify(expt.regions)
end

"""
    integrate!(region::Region, expt::Experiment1D)

Calculate integrals for a specific region across all slices.
"""
function integrate!(region::Region, expt::Experiment1D)
    @debug "Integrating region $(region.label[])"
    
    # Get the number of slices
    nslices = length(expt.specdata.y)
    
    # Allocate array for integrals
    integrals = zeros(nslices)
    
    # Iterate through slices
    for i in 1:nslices
        # Get x and y data for this slice
        x = expt.specdata.x[i]
        y = expt.specdata.y[i]
        
        # Find points within the region
        idx = findall(x -> region.xstart[] <= x <= region.xend[], x)
        
        # Calculate integral
        if !isempty(idx)
            # Calculate integral using trapezoidal rule
            dx = diff(x[idx])
            ymean = (y[idx[1:end-1]] + y[idx[2:end]]) ./ 2
            integrals[i] = sum(dx .* ymean)
        end
    end
    
    # Update the parameter value
    if haskey(region.parameters, :integral)
        if length(integrals) == 1
            region.parameters[:integral].value[] = integrals[1]
        else
            region.parameters[:integral].value[] = integrals
        end
    end
    
    # Mark as touched
    region.touched[] = true
    
    return integrals
end

"""
    fit!(expt::Experiment1D)

Fit all touched regions in the experiment.
"""
function fit!(expt::Experiment1D)
    # Check if anything has changed
    anythingchanged = any(region -> region.touched[], expt.regions[])
    
    if !anythingchanged || !expt.isfitting[]
        return
    end
    
    @async begin
        expt.state[][:mode][] = :fitting
        sleep(0.1) # Allow time for mode change to be processed
    end
    
    @async begin # Do fitting in a separate task
        sleep(0.1) # Allow time for mode change to be processed
        
        # Integrate all regions first
        integrate!(expt)
        
        # Fit each touched region
        for (i, region) in enumerate(expt.regions[])
            if region.touched[]
                fit!(region, expt)
                postfit!(region, expt)
            end
        end
        
        # Notify that regions have been updated
        notify(expt.regions)
        
        expt.state[][:mode][] = :normal
    end
end

"""
    fit!(region::Region, expt::Experiment1D)

Fit a specific region based on the experiment type.
"""
function fit!(region::Region, expt::Experiment1D)
    @debug "Fitting region $(region.label[])"
    
    # This is a placeholder - implementation will vary by experiment type
    # Each experiment type should override this method
    
    # Mark as no longer touched
    region.touched[] = false
end

"""
    setup_region_parameters!(region::Region, expt::Experiment1D)

Set up parameters for a region based on the experiment type.
"""
function setup_region_parameters!(region::Region, expt::Experiment1D)
    # Default implementation just adds the integral parameter
    region.parameters[:integral] = Parameter("Integral", 0.0)
end

"""
    ppm_range(region::Region)

Get the chemical shift range of a region.
"""
function ppm_range(region::Region)
    xstart, xend = minmax(region.xstart[], region.xend[])
    return xstart..xend
end

"""
    width(region::Region)

Get the width of a region in ppm.
"""
function width(region::Region)
    return abs(region.xend[] - region.xstart[])
end

"""
    center(region::Region)

Get the center point of a region.
"""
function center(region::Region)
    return (region.xstart[] + region.xend[]) / 2
end