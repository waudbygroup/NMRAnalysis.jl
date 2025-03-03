"""
    gui!(expt::Experiment1D)

Create and display the GUI for a 1D experiment.
"""
function gui!(expt::Experiment1D)
    GLMakie.activate!(; title="NMRAnalysis.jl 1D Analysis",
                      focus_on_show=true)

    state = expt.state[]
    
    g = Dict{Symbol,Any}()  # GUI state
    state[:gui][] = g
    
    # Create figure with layout
    g[:fig] = Figure(; size=(1200, 800))
    g[:paneltop] = g[:fig][1, 1:2] = GridLayout()
    g[:panelspectrum] = g[:fig][2:3, 1] = GridLayout()
    g[:panelinfo] = g[:fig][2, 2] = GridLayout()
    g[:panelplot] = g[:fig][3, 2] = GridLayout()
    
    # Top panel with controls
    create_top_panel!(g, state, expt)
    
    # Spectrum panel
    create_spectrum_panel!(g, state, expt)
    
    # Info panel
    create_info_panel!(g, state, expt)
    
    # Plot panel (experiment-specific)
    create_plot_panel!(g, state, expt)
    
    # Add event handlers
    add_event_handlers!(g, state, expt)
    
    # Display the figure
    display(g[:fig])
    
    return g[:fig]
end

"""
    create_top_panel!(g, state, expt)

Create the top panel with controls.
"""
function create_top_panel!(g, state, expt)
    # Buttons and controls
    g[:cmdresetzoom] = Button(g[:paneltop][1, 1]; label="Reset zoom")
    
    # Slice controls (only if more than one slice)
    if nslices(expt) > 1
        g[:sliderslice] = Slider(g[:paneltop][1, 2]; range=1:nslices(expt))
        g[:cmdsliceleft] = Button(g[:paneltop][1, 3]; label="←")
        g[:cmdsliceright] = Button(g[:paneltop][1, 4]; label="→")
        g[:slicelabel] = Label(g[:paneltop][1, 5], state[:current_slice_label])
    end
    
    # Fitting toggle
    g[:togglefit] = Toggle(g[:paneltop][1, 6]; active=true)
    Label(g[:paneltop][1, 7], "Fitting")
    
    # File operations
    g[:cmdload] = Button(g[:paneltop][1, 8]; label="Load regions")
    g[:cmdsave] = Button(g[:paneltop][1, 9]; label="Save results")
    
    # Custom controls for specific experiment types
    add_custom_controls!(g, state, expt)
end

"""
    create_spectrum_panel!(g, state, expt)

Create the main spectrum display panel.
"""
function create_spectrum_panel!(g, state, expt)
    # Create axis for spectrum
    g[:axspectrum] = Axis(g[:panelspectrum][1, 1];
                         xlabel="Chemical shift (ppm)",
                         ylabel="Intensity",
                         xreversed=true)
    
    # Plot the spectrum
    # Start with first slice
    xdata = expt.specdata.xplot
    ydata = expt.specdata.yplot
    
    # Plot spectrum
    g[:pltspectrum] = lines!(g[:axspectrum], xdata, ydata; color=:black)
    
    # Add visualization for integration regions
    g[:pltregions] = Dict{Int,Any}()
    
    # Add baseline
    g[:pltbaseline] = hlines!(g[:axspectrum], [0]; color=:gray, linestyle=:dash)
    
    # Auto-adjust limits initially
    on(g[:pltspectrum].input_args[2]) do _
        autolimits!(g[:axspectrum])
    end
end

"""
    create_info_panel!(g, state, expt)

Create the information panel.
"""
function create_info_panel!(g, state, expt)
    # Region operations
    g[:cmdrename] = Button(g[:panelinfo][1, 1]; label="(R)ename region")
    g[:cmddelete] = Button(g[:panelinfo][1, 2]; label="(D)elete region")
    
    # Help text
    Label(g[:panelinfo][2, 1:2], "Press (A) to add new region under cursor"; word_wrap=true)
    
    # Region info text
    g[:infotext] = Label(g[:panelinfo][3, 1:2], state[:current_region_info]; word_wrap=true)
    
    # Experiment info
    g[:exptinfo] = Label(g[:panelinfo][4, 1:2], experimentinfo(expt); word_wrap=true)
end

"""
    create_plot_panel!(g, state, expt)

Create the experiment-specific plot panel.
"""
function create_plot_panel!(g, state, expt)
    # Default implementation - creates a placeholder panel
    g[:axplot] = Axis(g[:panelplot][1, 1];
                     xlabel="x",
                     ylabel="y")
    
    # Specific implementation provided by experiment types
    create_experiment_plot!(g, state, expt)
end

"""
    create_experiment_plot!(g, state, expt)

Create the experiment-specific plot. Should be implemented by each experiment type.
"""
function create_experiment_plot!(g, state, expt)
    # Default implementation - empty plot
    # This should be overridden by specific experiment types
end

"""
    add_custom_controls!(g, state, expt)

Add experiment-specific controls. Should be implemented by each experiment type.
"""
function add_custom_controls!(g, state, expt)
    # Default implementation - no additional controls
    # This should be overridden by specific experiment types
end

"""
    add_event_handlers!(g, state, expt)

Set up event handlers for the GUI.
"""
function add_event_handlers!(g, state, expt)
    # Set up background color based on mode
    g[:fig].scene.backgroundcolor = lift(state[:mode]) do mode
        if mode == :fitting
            RGBAf(1.0, 0.63, 0.48, 1.0)      # :salmon
        elseif mode == :renaming || mode == :renamingstart
            RGBAf(0.75, 0.94, 1.0, 1.0)     # :lightblue
        elseif mode == :moving
            RGBAf(0.6, 0.98, 0.6, 1.0)      # :palegreen
        elseif mode == :resizing
            RGBAf(1.0, 0.87, 0.6, 1.0)      # :peachpuff
        else
            RGBAf(1.0, 1.0, 1.0, 1.0)       # :white
        end
    end
    
    # Reset zoom button
    on(g[:cmdresetzoom].clicks) do _
        reset_limits!(g[:axspectrum])
    end
    
    # Slice navigation if applicable
    if haskey(g, :sliderslice)
        connect!(state[:current_slice], g[:sliderslice].value)
        on(g[:cmdsliceleft].clicks) do _
            i = state[:current_slice][]
            if i > 1
                i -= 1
                set_close_to!(g[:sliderslice], i)
            end
        end
        on(g[:cmdsliceright].clicks) do _
            i = state[:current_slice][]
            if i < nslices(expt)
                i += 1
                set_close_to!(g[:sliderslice], i)
            end
        end
        
        # Update displayed spectrum when slice changes
        on(state[:current_slice]) do _
            update_spectrum_display!(expt)
        end
    end
    
    # Fitting toggle
    connect!(expt.isfitting, g[:togglefit].active)
    
    # Load regions
    on(g[:cmdload].clicks) do _
        loadregions!(expt)
    end
    
    # Save results
    on(g[:cmdsave].clicks) do _
        saveresults!(expt)
    end
    
    # Delete region
    on(g[:cmddelete].clicks) do _
        idx = state[:current_region_idx][]
        if idx > 0
            deleteregion!(expt, idx)
        end
    end
    
    # Rename region
    on(g[:cmdrename].clicks) do _
        if state[:current_region_idx][] > 0
            renameregion!(expt, state, :mouse)
        end
    end
    
    # Mouse handlers for regions
    on(events(g[:axspectrum]).mousebutton) do event
        process_mousebutton(expt, state, event)
    end
    
    on(events(g[:axspectrum]).mouseposition) do pos
        process_mouseposition(expt, state, pos)
    end
    
    # Keyboard handlers
    on(events(g[:axspectrum]).keyboardbutton) do event
        process_keyboardbutton(expt, state, event)
    end
    
    on(events(g[:fig]).unicode_input) do char
        process_unicode_input(expt, state, char)
    end
    
    # Add experiment-specific handlers
    add_experiment_handlers!(g, state, expt)
end

"""
    add_experiment_handlers!(g, state, expt)

Add experiment-specific event handlers. Should be implemented by each experiment type.
"""
function add_experiment_handlers!(g, state, expt)
    # Default implementation - no additional handlers
    # This should be overridden by specific experiment types
end

"""
    update_spectrum_display!(expt::Experiment1D)

Update the displayed spectrum data.
"""
function update_spectrum_display!(expt::Experiment1D)
    # Get current slice
    slice_idx = expt.state[][:current_slice][]
    
    # Update x and y data
    if slice_idx > 0 && slice_idx <= nslices(expt)
        expt.specdata.xplot[] = expt.specdata.x[slice_idx]
        expt.specdata.yplot[] = expt.specdata.y[slice_idx]
    end
end

"""
    update_region_display!(expt::Experiment1D)

Update the visualization of integration regions.
"""
function update_region_display!(expt::Experiment1D)
    g = expt.state[][:gui][]
    state = expt.state[]
    
    # Clear existing region visualizations
    for (_, plt) in g[:pltregions]
        delete!(plt)
    end
    empty!(g[:pltregions])
    
    # Add visualization for each region
    for (idx, region) in enumerate(expt.regions[])
        # Get region properties
        xstart = region.xstart[]
        xend = region.xend[]
        color = region.color[]
        
        # Get current slice y data
        slice_idx = state[:current_slice][]
        
        if slice_idx > 0 && slice_idx <= nslices(expt)
            x = expt.specdata.x[slice_idx]
            y = expt.specdata.y[slice_idx]
            
            # Find points within region
            idx_in_region = findall(xi -> xstart <= xi <= xend, x)
            
            if !isempty(idx_in_region)
                # Create region fill
                x_region = x[idx_in_region]
                y_region = y[idx_in_region]
                
                # Create polygon vertices for fill
                xs = [x_region; reverse(x_region)]
                ys = [y_region; zeros(length(y_region))]
                
                # Add region visualization
                g[:pltregions][idx] = poly!(g[:axspectrum], 
                                           Point2f.(xs, ys); 
                                           color=(color, 0.3), 
                                           strokewidth=2, 
                                           strokecolor=color)
                
                # Add vertical lines at region boundaries
                vlines!(g[:axspectrum], [xstart, xend]; 
                       color=color, 
                       linestyle=:dash)
                
                # Add label at top of region
                middle_x = (xstart + xend) / 2
                middle_y_idx = findfirst(xi -> xi >= middle_x, x_region)
                if !isnothing(middle_y_idx)
                    middle_y = maximum(y_region) * 1.05
                    text!(g[:axspectrum], 
                         middle_x, middle_y, 
                         text=region.label[], 
                         color=color,
                         align=(:center, :bottom))
                end
            end
        end
    end
end

"""
    create_experiment_specific_visualization!(g, state, expt, ::Type{<:Experiment1D})

Create experiment-specific visualizations based on the experiment type.
"""
function create_experiment_specific_visualization!(g, state, expt, ::Type{<:Experiment1D})
    # Default implementation - no special visualization
    # This should be overridden by specific experiment types
end