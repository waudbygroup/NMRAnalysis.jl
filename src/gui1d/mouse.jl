"""
    process_mousebutton(expt, state, event)

Process mouse button events.
"""
function process_mousebutton(expt, state, event)
    # Handle based on current mode
    if state[:mode][] == :normal && event.button == Mouse.left && event.action == Mouse.press
        # Left button pressed in normal mode
        handle_normal_click(expt, state, event)
    elseif (state[:mode][] == :moving || state[:mode][] == :resizing) && 
           event.button == Mouse.left && event.action == Mouse.release
        # Left button released while moving or resizing
        handle_end_drag(expt, state, event)
    end
    
    return Consume(false)
end

"""
    handle_normal_click(expt, state, event)

Handle mouse click in normal mode.
"""
function handle_normal_click(expt, state, event)
    g = state[:gui][]
    pos = events(g[:axspectrum]).mouseposition[]
    screen_pos = Point2f(pos)
    
    # Convert to data coordinates
    data_pos = data_coordinates(g[:axspectrum], screen_pos)
    x_pos = data_pos[1]
    y_pos = data_pos[2]
    
    # Store start position for potential dragging
    state[:mouse_start_pos][] = (x_pos, y_pos)
    
    # Check if clicked on a region edge
    edge, region_idx = find_region_edge(expt, state, x_pos)
    
    if edge != :none
        # Clicked on a region edge - start resizing
        state[:mode][] = :resizing
        state[:drag_region_idx][] = region_idx
        state[:drag_edge][] = edge
        state[:current_region_idx][] = region_idx
        return Consume(true)
    end
    
    # Check if clicked inside a region
    region_idx = find_region_at_position(expt, state, x_pos)
    
    if region_idx > 0
        # Clicked on a region - start moving
        state[:mode][] = :moving
        state[:drag_region_idx][] = region_idx
        state[:current_region_idx][] = region_idx
        return Consume(true)
    end
    
    # If 'a' key is pressed, add a new region
    if ispressed(g[:fig], Keyboard.a)
        addregion!(expt, x_pos)
        return Consume(true)
    end
    
    # No special action, just update current region selection
    state[:current_region_idx][] = 0
    return Consume(false)
end

"""
    handle_end_drag(expt, state, event)

Handle mouse release after dragging.
"""
function handle_end_drag(expt, state, event)
    # Reset mode and drag info
    state[:mode][] = :normal
    region_idx = state[:drag_region_idx][]
    state[:drag_region_idx][] = 0
    state[:drag_edge][] = :none
    
    # Update integration and fitting
    if region_idx > 0 && expt.isfitting[]
        # Mark region as touched to trigger fitting
        expt.regions[][region_idx].touched[] = true
        notify(expt.regions)
    end
    
    return Consume(true)
end

"""
    process_mouseposition(expt, state, mousepos)

Process mouse movement events.
"""
function process_mouseposition(expt, state, mousepos)
    g = state[:gui][]
    screen_pos = Point2f(mousepos)
    
    # Convert to data coordinates
    data_pos = data_coordinates(g[:axspectrum], screen_pos)
    x_pos = data_pos[1]
    y_pos = data_pos[2]
    
    # Update current position
    state[:mouse_current_pos][] = (x_pos, y_pos)
    
    # Handle based on current mode
    if state[:mode][] == :moving
        # Moving a region
        handle_region_move(expt, state)
        return Consume(true)
    elseif state[:mode][] == :resizing
        # Resizing a region
        handle_region_resize(expt, state)
        return Consume(true)
    elseif state[:mode][] == :normal
        # Normal mode - update cursor based on what's under it
        update_cursor(expt, state, x_pos)
    end
    
    return Consume(false)
end

"""
    handle_region_move(expt, state)

Handle dragging to move a region.
"""
function handle_region_move(expt, state)
    region_idx = state[:drag_region_idx][]
    if region_idx <= 0 || region_idx > nregions(expt)
        return
    end
    
    # Calculate delta x from start of drag
    start_x = state[:mouse_start_pos][][1]
    current_x = state[:mouse_current_pos][][1]
    delta_x = current_x - start_x
    
    # Get region
    region = expt.regions[][region_idx]
    
    # Update region position
    old_start = region.xstart[]
    old_end = region.xend[]
    region.xstart[] = old_start + delta_x
    region.xend[] = old_end + delta_x
    
    # Update start position for next movement
    state[:mouse_start_pos][] = state[:mouse_current_pos][]
    
    # Mark as touched and notify
    region.touched[] = true
    notify(expt.regions)
    
    # Update display
    update_region_display!(expt)
end

"""
    handle_region_resize(expt, state)

Handle dragging to resize a region.
"""
function handle_region_resize(expt, state)
    region_idx = state[:drag_region_idx][]
    if region_idx <= 0 || region_idx > nregions(expt)
        return
    end
    
    # Get current x position
    current_x = state[:mouse_current_pos][][1]
    
    # Get region
    region = expt.regions[][region_idx]
    
    # Update region based on which edge is being dragged
    edge = state[:drag_edge][]
    if edge == :left
        region.xstart[] = current_x
    elseif edge == :right
        region.xend[] = current_x
    end
    
    # Ensure xstart < xend
    if region.xstart[] > region.xend[]
        xstart = region.xend[]
        xend = region.xstart[]
        region.xstart[] = xstart
        region.xend[] = xend
        
        # Flip edge being dragged
        if edge == :left
            state[:drag_edge][] = :right
        elseif edge == :right
            state[:drag_edge][] = :left
        end
    end
    
    # Mark as touched and notify
    region.touched[] = true
    notify(expt.regions)
    
    # Update display
    update_region_display!(expt)
end

"""
    find_region_edge(expt, state, x_pos)

Find if the given x position is near a region edge.
Returns (edge_type, region_index) where edge_type is :none, :left, or :right.
"""
function find_region_edge(expt, state, x_pos)
    # Define edge sensitivity in data units (ppm)
    edge_sensitivity = 0.05
    
    for (idx, region) in enumerate(expt.regions[])
        # Check left edge
        if abs(region.xstart[] - x_pos) <= edge_sensitivity
            return :left, idx
        end
        
        # Check right edge
        if abs(region.xend[] - x_pos) <= edge_sensitivity
            return :right, idx
        end
    end
    
    return :none, 0
end

"""
    find_region_at_position(expt, state, x_pos)

Find which region contains the given x position.
Returns region index or 0 if none.
"""
function find_region_at_position(expt, state, x_pos)
    for (idx, region) in enumerate(expt.regions[])
        if region.xstart[] <= x_pos <= region.xend[]
            return idx
        end
    end
    
    return 0
end

"""
    update_cursor(expt, state, x_pos)

Update the mouse cursor based on what's under it.
"""
function update_cursor(expt, state, x_pos)
    g = state[:gui][]
    
    # Check if over a region edge
    edge, _ = find_region_edge(expt, state, x_pos)
    
    if edge == :left || edge == :right
        # Show resize cursor
        GLMakie.set_cursor(g[:fig].scene, :ew_resize) # East-west resize cursor
    else
        # Check if inside a region
        region_idx = find_region_at_position(expt, state, x_pos)
        
        if region_idx > 0
            # Show move cursor
            GLMakie.set_cursor(g[:fig].scene, :move)
        else
            # Default cursor
            GLMakie.set_cursor(g[:fig].scene, :arrow)
        end
    end
end

"""
    data_coordinates(ax, screen_pos)

Convert screen coordinates to data coordinates.
"""
function data_coordinates(ax, screen_pos)
    # Get transformation from screen to data coordinates
    transform = Makie.inv_transform(ax)
    
    # Apply transformation
    data_pos = transform(screen_pos)
    
    return data_pos
end