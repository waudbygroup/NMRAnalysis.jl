"""
    process_keyboardbutton(expt, state, event)

Process keyboard button events.
"""
function process_keyboardbutton(expt, state, event)
    @debug "Keyboard event: $event"
    g = state[:gui][]
    
    if state[:mode][] == :normal && event.action == Keyboard.press
        # Normal mode keyboard shortcuts
        handle_normal_keypress(expt, state, event)
    elseif (state[:mode][] == :renaming || state[:mode][] == :renamingstart) && event.action == Keyboard.press
        # Renaming mode keyboard handling
        handle_renaming_keypress(expt, state, event)
    end
end

"""
    handle_normal_keypress(expt, state, event)

Handle keyboard press in normal mode.
"""
function handle_normal_keypress(expt, state, event)
    g = state[:gui][]
    
    if ispressed(g[:fig], Keyboard.a)
        # Add new region at mouse position
        x_pos = state[:mouse_current_pos][][1]
        addregion!(expt, x_pos)
    elseif ispressed(g[:fig], Keyboard.d)
        # Delete current region
        idx = state[:current_region_idx][]
        if idx > 0
            state[:current_region_idx][] = 0
            deleteregion!(expt, idx)
        end
    elseif ispressed(g[:fig], Keyboard.r)
        # Rename current region
        if state[:current_region_idx][] > 0
            renameregion!(expt, state, :keyboard)
        end
    elseif haskey(g, :sliderslice) && ispressed(g[:fig], Keyboard.left)
        # Previous slice
        i = state[:current_slice][]
        if i > 1
            i -= 1
            set_close_to!(g[:sliderslice], i)
        end
    elseif haskey(g, :sliderslice) && ispressed(g[:fig], Keyboard.right)
        # Next slice
        i = state[:current_slice][]
        if i < nslices(expt)
            i += 1
            set_close_to!(g[:sliderslice], i)
        end
    end
end

"""
    handle_renaming_keypress(expt, state, event)

Handle keyboard press in renaming mode.
"""
function handle_renaming_keypress(expt, state, event)
    if event.key == Keyboard.enter
        # Complete renaming
        region = state[:current_region][]
        region.label[] = region.label[][1:(end - 1)]  # Remove cursor character
        state[:mode][] = :normal
        notify(expt.regions)
        return Consume()
    elseif event.key == Keyboard.backspace
        # Delete character
        region = state[:current_region][]
        if length(region.label[]) > 1
            region.label[] = region.label[][1:(end - 2)] * "‸"
            notify(expt.regions)
        end
        return Consume()
    elseif event.key == Keyboard.escape
        # Cancel renaming
        region = state[:current_region][]
        region.label[] = state[:oldlabel][]
        notify(expt.regions)
        state[:mode][] = :normal
        return Consume()
    end
    
    return Consume(false)
end

"""
    process_unicode_input(expt, state, character)

Process unicode input (for text entry).
"""
function process_unicode_input(expt, state, character)
    @debug "Processing unicode input: $character"
    
    if state[:mode][] == :renamingstart
        state[:mode][] = :renaming
        if character == 'r'
            # Discard initial 'r' character from keyboard shortcut
            return Consume()
        end
    end
    
    if state[:mode][] == :renaming
        # Add character to label
        region = state[:current_region][]
        region.label[] = region.label[][1:(end - 1)] * character * "‸"
        notify(expt.regions)
        return Consume()
    end
    
    return Consume(false)
end