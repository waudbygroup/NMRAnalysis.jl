"""
    preparestate(expt::Experiment)

Initialize the state dictionary for an experiment.
"""
function preparestate!(expt::Experiment)
    @debug "Preparing state"
    state = expt.state

    # Mode can be :normal, :renaming, :renamingstart, :moving, :resizing, :fitting
    state[:mode] = Observable(:normal)

    # Counter for total regions created
    state[:total_regions] = Observable(0)

    # Current slice index
    state[:current_slice] = Observable(1)
    
    # Current slice label
    state[:current_slice_label] = lift(idx -> slicelabel(expt, idx), state[:current_slice])

    # Current region index
    state[:current_region_idx] = Observable(0)
    
    # Current region
    state[:current_region] = Observable{Union{Region,Nothing}}(nothing)
    
    # Map current_region based on current_region_idx
    map!(state[:current_region], expt.regions, state[:current_region_idx]) do regions, idx
        if idx > 0 && idx <= length(regions)
            regions[idx]
        else
            nothing
        end
    end
    
    # Region info text
    state[:current_region_info] = lift(idx -> regioninfo(expt, idx), state[:current_region_idx])
    
    # When renaming, store the old label
    state[:oldlabel] = Observable("")
    
    # Initialize GUI reference
    state[:gui] = Observable{Dict{Symbol,Any}}(Dict{Symbol,Any}())
    
    # Add experiment-specific state
    completestate!(state, expt)
    
    return expt
end

"""
    completestate!(state, expt::Experiment)

Complete initialization of state with experiment-specific observables.
"""
function completestate!(state, ::Experiment)
    # Default implementation - no additional state needed
    return state
end