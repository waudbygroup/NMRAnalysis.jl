function process_keyboardbutton(state, event)
    if state[:renaming][]
        if event.action == Keyboard.press && event.key == Keyboard.enter
            state[:peaks][][state[:activepeakindex]].label = state[:peaks][][state[:activepeakindex]].label[1:end-1]
            state[:renaming][] = false
            notify(state[:peaks])
            return Consume()
        elseif event.action == Keyboard.press && event.key == Keyboard.backspace
            if length(state[:peaks][][state[:activepeakindex]].label) > 1
                state[:peaks][][state[:activepeakindex]].label = state[:peaks][][state[:activepeakindex]].label[1:end-2] * "‸"
                notify(state[:peaks])
                return Consume()
            end
        elseif event.action == Keyboard.press && event.key == Keyboard.escape
            # restore previous label
            state[:peaks][][state[:activepeakindex]].label = state[:oldlabel]
            notify(state[:peaks])

            state[:renaming][] = false
            return(Consume())
        end
        return Consume(false)
    end

    # not renaming - parse for peak commands
    event.action == Keyboard.press || return Consume(false)

    if event.key == Keyboard.a
        # Add marker
        newpos = mouseposition(state[:gui][:mainax].scene)
        
        newpeak = createpeak!(state, newpos)
        addpeakandrecluster!(state, newpeak)
        return Consume()
    end

    pick_plt, pick_i = pick(state[:gui][:fig].scene, Makie.mouseposition_px(state[:gui][:fig].scene), 50)
    pick_plt == state[:gui][:peakplot] || return Consume(false)

    if event.key == Keyboard.d
        # Delete marker
        deletepeakandrecluster!(state, pick_i)
        return Consume()
    elseif event.key == Keyboard.r
        # Rename marker
        state[:activepeakindex] = pick_i
        state[:oldlabel] = state[:peaks][][state[:activepeakindex]].label # store previous label
        state[:peaks][][state[:activepeakindex]].label = "‸" # clear label
        state[:renaming][] = true
        state[:startinglabelinput] = true

        notify(state[:peaks])
        return Consume()
    end

    return Consume(false)
end


function process_unicode_input(state, character)
    if state[:renaming][]
        if state[:startinglabelinput] && character == 'r'
            # discard 'r' character hanging over from initial keypress
            state[:startinglabelinput] = false
            return Consume()
        end
        state[:peaks][][state[:activepeakindex]].label = state[:peaks][][state[:activepeakindex]].label[1:end-1] * character * "‸"
        notify(state[:peaks])
        return Consume()
    end
    return Consume(false)
end