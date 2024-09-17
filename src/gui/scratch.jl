using Revise

using DataStructures
using GLMakie
using LightGraphs
using LsqFit
using NMRTools

include("util.jl")
include("types.jl")
include("models.jl")
include("peaks.jl")
include("specdata.jl")
include("keyboard.jl")
include("mouse.jl")
include("mask.jl")
include("clustering.jl")
include("sim.jl")
include("fitting.jl")

spectra = [
    loadnmr("/Users/chris/NMR/crick-950/pallavi_trxl2_240322/4/"),
    # loadnmr("/Users/chris/NMR/crick-800/jiwoo_CTD_his_20240725/804"),
    # loadnmr("/Users/chris/NMR/crick-800/jiwoo_CTD_his_20240725/806"),
]
model = ModelAmplitudes()


# create dictionaries to hold the data, state, and gui elements
state = Dict()
state[:gui] = Dict()
state[:specdata] = preparespecdata(spectra, 6 .. 10)
gui = state[:gui] # for convenience
specdata = state[:specdata] # for convenience



# now create a state dictionary
state[:model] = model
state[:slice] = Observable(1)
state[:x] = Observable(specdata[:x][1])
state[:y] = Observable(specdata[:y][1])
state[:z] = Observable(specdata[:z][1])
state[:zfit] = Observable(specdata[:zfit][1])
state[:mask] = Observable(nothing)
state[:clev0] = Observable(5.0)
state[:cfactor] = 1.7
state[:baselev] = vcat(-1 * state[:cfactor] .^ (6:-1:0), state[:cfactor] .^ (0:6))
state[:peaks] = Observable(Vector{Peak2D}())
state[:peakcounter] = 0
state[:renaming] = Observable(false)
state[:dragging] = Observable(false)
state[:highlighted] = Observable(0)
state[:Xradius] = Observable(0.06)
state[:Yradius] = Observable(0.04)
state[:fitpositions] = Observable(Vector{Point2f}())
state[:fixpeakpositions] = fixpeakpositions(model)
state[:fixamplitudes] = fixamplitudes(model)
state[:fixlinewidths] = fixlinewidths(model)

# create the figure
gui[:fig] = Figure()
gui[:mainpanel] = gui[:fig][1, 1]
gui[:sidepanel] = gui[:fig][1, 2]

# create the control panel on the right hand side
gui[:sg_spectrum] = SliderGrid(
    gui[:sidepanel][1, 1],
    (label = "Spectrum", range = 1:specdata[:nspec], format = "{}", startvalue = 1),
    tellheight = false)
gui[:slider_spectrum] = gui[:sg_spectrum].sliders[1]
connect!(state[:slice], gui[:slider_spectrum].value)
gui[:cmd_lower] = Button(gui[:sidepanel][2,1][1,1], label="Lower Levels")
gui[:cmd_raise] = Button(gui[:sidepanel][2,1][1,2], label="Raise Levels")

gui[:sg_radii] = SliderGrid(
    gui[:sidepanel][3, 1],
    (label = "X fitting radius", range = 0.005:0.005:0.1, format = "{:.3f}", startvalue = 0.06),
    (label = "Y fitting radius", range = 0.05:0.05:1, format = "{:.2f}", startvalue = 0.4),
    tellheight = false)
gui[:slider_Xradius] = gui[:sg_radii].sliders[1]
gui[:slider_Yradius] = gui[:sg_radii].sliders[2]
connect!(state[:Xradius], gui[:slider_Xradius].value)
connect!(state[:Yradius], gui[:slider_Yradius].value)

gui[:toggle_fit] = Toggle(gui[:sidepanel][4,1][1,1], active=false)
Label(gui[:sidepanel][4,1][1,2], "Fitting")

gui[:cmd_quit] = Button(gui[:sidepanel][5,1], label="Quit")

# callbacks
on(state[:slice]) do i
    # update the axis values and data simultaneously:
    state[:x].val = specdata[:x][i]
    state[:y].val = specdata[:y][i]
    state[:zfit][] = specdata[:zfit][i]
    state[:z][] = specdata[:z][i]
    state[:fitpositions][] = getfitpositions(state)
end
on(gui[:cmd_lower].clicks) do i
    state[:clev0][] /= state[:cfactor]
end
on(gui[:cmd_raise].clicks) do i
    state[:clev0][] *= state[:cfactor]
end
on(gui[:cmd_quit].clicks) do i
    GLMakie.GLFW.SetWindowShouldClose(gui[:window], true) # this will close the window after all callbacks are finished
end

# now create the contour plot panel
gui[:mainax] = Axis(gui[:mainpanel],
    xreversed=true, yreversed=true,
    xgridvisible=false, ygridvisible=false,
    xlabel=@lift(specdata[:xlabels][$(state[:slice])]),
    ylabel=@lift(specdata[:ylabels][$(state[:slice])]),
    title=@lift(specdata[:titles][$(state[:slice])])
)
maskplot = heatmap!(gui[:mainax], state[:x], state[:y], getmask(state),
    colormap=[colorant"white", colorant"navajowhite"], colorrange=(0,1), inspectable=false)
gui[:contourplot] = contour!(gui[:mainax], state[:x], state[:y], state[:z],
    levels=@lift(state[:baselev] * $(state[:clev0])),
    colormap=[colorant"pink", colorant"grey"],
    colorrange=(-0.001, 0.001),
    lowclip=:pink,
    highclip=:grey,
    )
gui[:fitplot] = contour!(gui[:mainax], state[:x], state[:y], state[:zfit],
    levels=@lift(state[:baselev] * $(state[:clev0])),
    colormap=[colorant"blue", colorant"red"],
    colorrange=(-0.001, 0.001),
    lowclip=:blue,
    highclip=:red,
    )
gui[:fitpeakplot] = scatter!(gui[:mainax], state[:fitpositions],
    color=:darkblue,
    markersize=15,
    marker=:+,
    # inspectable=true,
    )
gui[:peakplot] = scatter!(gui[:mainax], getpeakpositions(state),
        color=getpeakcolors(state),
        markersize=15,
        marker=:x,
        # inspectable=true,
        )
gui[:labelsplot] = text!(gui[:mainax], getpeaklabels(state), offset=(10,10))
connect!(gui[:fitplot].visible, gui[:toggle_fit].active)
connect!(gui[:fitpeakplot].visible, gui[:toggle_fit].active)


# key commands
on(events(gui[:mainax]).keyboardbutton, priority = 2) do event
    process_keyboardbutton(state, event)
end
on(events(gui[:fig]).unicode_input) do character
    process_unicode_input(state, character)
end

# mouse dragging
on(events(gui[:mainax]).mousebutton, priority = 2) do event
    process_mousebutton(state, event)
end
on(events(gui[:mainax]).mouseposition, priority = 2) do mousepos
    process_mouseposition(state, mousepos)
end

# fitting
onany(state[:peaks], gui[:toggle_fit].active) do _, fitting
    fitting || return Consume(false)
    state[:dragging][] && return Consume(false)

    fit!(state)
    return Consume()
end

# connect background colour to state
bgcolor = lift(x -> x ? colorant"lightgreen" : colorant"white", state[:renaming])
connect!(gui[:fig].scene.backgroundcolor, bgcolor)

# display the window
gui[:window] = GLMakie.to_native(display(gui[:fig]))


