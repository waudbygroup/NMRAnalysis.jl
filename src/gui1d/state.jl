"""
    preparestate(expt::Experiment1D)

Initialize the state dictionary for an experiment.
"""
function preparestate(expt::Experiment1D)
    @debug "Preparing state"
    state = Dict{Symbol,Observable}()

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
    state[:current_region_info] = lift(idx -> regioninfotext(expt, idx), state[:current_region_idx])
    
    # When renaming, store the old label
    state[:oldlabel] = Observable("")
    
    # Information about plot limits
    state[:x_lims] = Observable((0.0, 10.0))
    state[:y_lims] = Observable((0.0, 1.0))
    
    # Mouse position for drag operations
    state[:mouse_start_pos] = Observable((0.0, 0.0))
    state[:mouse_current_pos] = Observable((0.0, 0.0))
    
    # Drag operation details
    state[:drag_region_idx] = Observable(0)
    state[:drag_edge] = Observable(:none)  # :none, :left, :right, :both
    
    # Initialize GUI reference
    state[:gui] = Observable{Dict{Symbol,Any}}(Dict{Symbol,Any}())
    
    # Add experiment-specific state
    completestate!(state, expt)
    
    return state
end

"""
    completestate!(state, expt::Experiment1D)

Complete initialization of state with experiment-specific observables.
"""
function completestate!(state, expt::Experiment1D)
    # Default implementation - no additional state needed
    return state
end