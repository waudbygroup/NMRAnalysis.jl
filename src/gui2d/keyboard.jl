function process_keyboardbutton(expt, state, event)
    @debug "keyboard event: $event"
    g = state[:gui][]
    if state[:mode][] == :normal
        if event.action == Keyboard.press && ispressed(g[:fig], Keyboard.a)
            pos = mouseposition(g[:axcontour])
            state[:total_peaks][] += 1
            id = "X$(state[:total_peaks][])"
            addpeak!(expt, Point2f(pos), id)
        elseif event.action == Keyboard.press && ispressed(g[:fig], Keyboard.d)
            idx = state[:current_peak_idx][]
            if idx > 0
                state[:current_peak_idx][] = 0
                deletepeak!(expt, idx)
            end
        elseif event.action == Keyboard.press && ispressed(g[:fig], Keyboard.r)
            if state[:current_peak_idx][] > 0
                renamepeak!(expt, state, :keyboard)
            end
        end
    elseif state[:mode][] == :renaming || state[:mode][] == :renamingstart
        if event.action == Keyboard.press && event.key == Keyboard.enter
            state[:current_peak][].label[] = state[:current_peak][].label[][1:end-1]
            state[:mode][] = :normal
            notify(expt.peaks)
            return Consume()
        elseif event.action == Keyboard.press && event.key == Keyboard.backspace
            if length(state[:peaks][][state[:activepeakindex]].label) > 1
                state[:current_peak][].label[] = state[:current_peak][].label[][1:end-2] * "‸"
                notify(expt.peaks)
                return Consume()
            end
        elseif event.action == Keyboard.press && event.key == Keyboard.escape
            # restore previous label
            state[:current_peak][].label[] = state[:oldlabel][]
            notify(expt.peaks)

            state[:mode][] = :normal
            return(Consume())
        end
        return Consume(false)
    end
end

function process_unicode_input(expt, state, character)
    @debug "Processing unicode input: $character"
    if state[:mode][] == :renamingstart
        state[:mode][] = :renaming
        if character == 'r'
            # discard 'r' character hanging over from initial keypress
            return Consume()
        end
    end
    if state[:mode][] == :renaming
        state[:current_peak][].label[] = state[:current_peak][].label[][1:end-1] * character * "‸"
        notify(expt.peaks)
        return Consume()
    end
    return Consume(false)
end