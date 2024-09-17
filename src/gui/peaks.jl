mutable struct Peak2D
    initial_position::MaybeVector{Point2f}
    # amplitude::MaybeVector{Float64}
    # lwX::MaybeVector{Float64}
    # lwY::MaybeVector{Float64}

    label::String
    touched::Bool

    pars#::Dict{Symbol, Float64}
    # pars_err::Dict{Symbol, Float64}
    sim_parameters
end
Peak2D(position, label="") = Peak2D(MaybeVector(position), label, true, DefaultDict{Symbol, Any}(nothing), DefaultDict{Symbol, Any}(nothing))

getpeakpositions(state) = @lift([p.initial_position[$(state[:slice])] for p ∈ $(state[:peaks])])
function getfitpositions(state)
    peaks = state[:peaks][]
    slice = state[:slice][]
    positions = map(peaks) do peak
        isnothing(peak.sim_parameters[:position]) ? nothing : peak.sim_parameters[:position][slice]
    end
    if any(isnothing.(positions)) || length(positions) == 0
        Vector{Point2f}()
    else
        positions
    end
end

function getpeakcolors(state)
    lift(state[:peaks], state[:highlighted]) do P, idx
        cols = [p.touched ? "orangered" : "dodgerblue" for p ∈ P]
        if idx > 0 && idx <= length(P)
            cols[idx] = "limegreen"
        end
        cols
    end
end

function getpeaklabels(state)
    lift(state[:peaks], state[:slice]) do P, slice
        if isempty(P)
            [("",
             Point2f((state[:x][][1]+state[:x][][end])/2,
                (state[:y][][1]+state[:y][][end])/2))]
        else
            [(p.label, p.initial_position[slice]) for p ∈ P]
        end
    end
end


function createpeak!(state, newpos)
    state[:peakcounter] += 1
    newlabel = "X$(state[:peakcounter])"
    positionpoint = Point2f(newpos)
    createpeak!(state, positionpoint, newlabel, state[:model])
    # if state[:fixpeakpositions]
    #     Peak2D(Point2f(newpos), newlabel)
    # else
    #     n = state[:specdata][:nspec]
    #     Peak2D([Point2f(newpos) for i=1:n], newlabel)
    # end
end

function addpeakandrecluster!(state, newpeak)
    peaks = state[:peaks]
    push!(peaks.val, newpeak)
    i = length(peaks[]) # index of new peak
    state[:clusters] = makeclusters(state)
    touchpeak!(peaks, i, state[:clusters])

    state[:highlighted][] = i # highlight new peak
    notify(state[:peaks])
    return i
end

function deletepeakandrecluster!(state, peakindex)
    oldclusters = state[:clusters]
    peaks = state[:peaks]
    touchpeak!(peaks, peakindex, oldclusters)
    deleteat!(peaks.val, peakindex)
    state[:clusters] = makeclusters(state)
    
    state[:highlighted][] = 0 # remove highlighting
    notify(state[:peaks])
end

function updatepeakandrecluster!(state, peakindex)
    oldclusters = state[:clusters]
    peaks = state[:peaks]
    touchpeak!(peaks, peakindex, oldclusters)
    state[:clusters] = makeclusters(state)
    touchpeak!(peaks, peakindex, state[:clusters])

    notify(state[:peaks])
end

function touchpeak!(peaks, peakindex, clusters)
    for cluster ∈ clusters
        if peakindex ∈ cluster
            for i ∈ cluster
                peaks[][i].touched = true
            end
            break
        end
    end
end

