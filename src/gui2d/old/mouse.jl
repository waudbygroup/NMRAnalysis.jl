function process_mousebutton(state, event)
    # consume input if in process of renaming a peak
    state[:renaming][] && return Consume()

    if event.button == Mouse.left && event.action == Mouse.press
        pick_plt, pick_i = pick(state[:gui][:fig].scene, Makie.mouseposition_px(state[:gui][:fig].scene), 50)
        state[:dragging][] = (pick_plt == state[:gui][:peakplot])
        if state[:dragging][]
            @debug "event: initiate drag" pick_i
            state[:dragpeakindex] = pick_i
        end
        return Consume(state[:dragging][])
    elseif event.button == Mouse.left && event.action == Mouse.release
        if state[:dragging][]
            @debug "event: stop dragging"
            state[:dragging][] = false
            # TODO
            updatepeakandrecluster!(state, state[:dragpeakindex])
            @debug "event: stop dragging - done"
        end
        return Consume(false)
    end
    return Consume(false)
end



function process_mouseposition(state, mousepos)
    state[:renaming][] && return Consume(true)

    if state[:dragging][]
        newpos = mouseposition(state[:gui][:mainax].scene)
        state[:peaks][][state[:dragpeakindex]].initial_position[state[:slice][]] = Point2f(newpos)
        state[:peaks][][state[:dragpeakindex]].touched = true

        notify(state[:peaks])
        return Consume()
    else
        pick_plt, pick_i = pick(state[:gui][:fig].scene, Makie.mouseposition_px(state[:gui][:fig].scene), 50)
        if(pick_plt == state[:gui][:peakplot])
            state[:highlighted][] = pick_i
            updatefitplot(pick_i, state)
        else
            state[:highlighted][] = 0
        end
    end
    return Consume(false)
end