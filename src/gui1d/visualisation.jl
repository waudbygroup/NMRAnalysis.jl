"""
    visualisationtype(expt::Experiment1D)

Get the visualization strategy for an experiment.
"""
visualisationtype(::Experiment1D) = LinePlotVisualisation()

"""
    create_experiment_plot!(g, state, expt::Experiment1D)

Create the experiment-specific plot based on visualization type.
"""
function create_experiment_plot!(g, state, expt::Experiment1D)
    create_experiment_plot!(g, state, expt, visualisationtype(expt))
end

"""
    create_experiment_plot!(g, state, expt, ::LinePlotVisualisation)

Create a line plot visualization.
"""
function create_experiment_plot!(g, state, expt, ::LinePlotVisualisation)
    g[:axplot] = Axis(g[:panelplot][1, 1];
                     xlabel="x",
                     ylabel="y")
    
    # Add placeholder for data
    g[:pltdata] = Observable(Point2f[])
    g[:plterrors] = Observable(Vector{Tuple{Float64,Float64,Float64}}[])
    g[:pltfit] = Observable(Point2f[])
    
    # Add plot elements
    scatterlines!(g[:axplot], g[:pltdata])
    errorbars!(g[:axplot], g[:plterrors]; whiskerwidth=10)
    lines!(g[:axplot], g[:pltfit]; color=:red, linewidth=2)
    
    # Set up update callback
    on(state[:current_region]) do region
        update_plot_data!(g, state, expt, region)
    end
end

"""
    create_experiment_plot!(g, state, expt, ::StackedPlotVisualisation)

Create a stacked plot visualization.
"""
function create_experiment_plot!(g, state, expt, ::StackedPlotVisualisation)
    g[:axplot] = Axis(g[:panelplot][1, 1];
                     xlabel="x",
                     ylabel="y")
    
    # Specific implementation for stacked plots
    # Similar to line plot but with vertical offsets
end

"""
    create_experiment_plot!(g, state, expt, ::WaterfallPlotVisualisation)

Create a 3D waterfall plot visualization.
"""
function create_experiment_plot!(g, state, expt, ::WaterfallPlotVisualisation)
    g[:axplot] = Axis3(g[:panelplot][1, 1];
                      xlabel="x",
                      ylabel="Slice",
                      zlabel="y")
    
    # Specific implementation for 3D waterfall plots
end

"""
    update_plot_data!(g, state, expt, region)

Update the plot data based on the selected region.
"""
function update_plot_data!(g, state, expt, region)
    if isnothing(region)
        # No region selected - clear plot
        g[:pltdata][] = Point2f[]
        g[:plterrors][] = Tuple{Float64,Float64,Float64}[]
        g[:pltfit][] = Point2f[]
        return
    end
    
    # Get model data from experiment
    update_model_visualization!(g, state, expt, region)
end

"""
    update_model_visualization!(g, state, expt, region)

Update the model visualization based on the experiment type.
This should be implemented by each experiment type.
"""
function update_model_visualization!(g, state, expt, region)
    # Default implementation - empty plot
    g[:pltdata][] = Point2f[]
    g[:plterrors][] = Tuple{Float64,Float64,Float64}[]
    g[:pltfit][] = Point2f[]
end

"""
    update_diffusion_visualization!(g, state, expt, region)

Update visualization for diffusion experiment.
"""
function update_model_visualization!(g, state, expt::DiffusionExperiment, region)
    # Get diffusion model data
    obs_points, obs_errors, fit_points = get_model_data(region, expt, DiffusionModel(expt.γ, expt.δ, expt.Δ, expt.σ))
    
    # Update plot data
    g[:pltdata][] = obs_points
    g[:plterrors][] = obs_errors
    g[:pltfit][] = fit_points
    
    # Update axis labels
    g[:axplot].xlabel = "Gradient strength / T m⁻¹"
    g[:axplot].ylabel = "Normalized integral"
end

"""
    update_tract_visualization!(g, state, expt, region)

Update visualization for TRACT experiment.
"""
function update_model_visualization!(g, state, expt::TRACTExperiment, region)
    # Get TRACT model data
    tract_data = get_tract_data(region, expt)
    if isempty(tract_data[1])
        # No data
        g[:pltdata][] = Point2f[]
        g[:plterrors][] = Tuple{Float64,Float64,Float64}[]
        g[:pltfit][] = Point2f[]
        return
    end
    
    # Extract data
    trosy_points, trosy_errors, trosy_fit, 
    antitrosy_points, antitrosy_errors, antitrosy_fit = tract_data
    
    # Combine TROSY and anti-TROSY data
    all_points = vcat(trosy_points, antitrosy_points)
    all_errors = vcat(trosy_errors, antitrosy_errors)
    all_fits = vcat(trosy_fit, antitrosy_fit)
    
    # Update plot data
    g[:pltdata][] = all_points
    g[:plterrors][] = all_errors
    g[:pltfit][] = all_fits
    
    # Update axis labels
    g[:axplot].xlabel = "Relaxation time / s"
    g[:axplot].ylabel = "Normalized integral"
    
    # Add legend
    if !haskey(g, :pltlegend)
        g[:pltlegend] = Legend(g[:panelplot][1, 2],
                             ["TROSY", "Anti-TROSY"],
                             [PolyElement(color=:blue), PolyElement(color=:red)])
    end
end

"""
    save_results_plots!(expt::Experiment1D, folder::AbstractString)

Save plots for all regions to the specified folder.
"""
function save_results_plots!(expt::Experiment1D, folder::AbstractString)
    CairoMakie.activate!()
    
    # Save overview plot
    save_overview_plot!(expt, folder)
    
    # Save individual region plots
    for region in expt.regions[]
        save_region_plot!(expt, region, folder)
    end
    
    # Switch back to interactive mode
    GLMakie.activate!()
end

"""
    save_overview_plot!(expt::Experiment1D, folder::AbstractString)

Save an overview plot of all regions.
"""
function save_overview_plot!(expt::Experiment1D, folder::AbstractString)
    fig = Figure(; size=(1200, 800))
    
    # Create spectrum plot
    ax = Axis(fig[1, 1];
             xlabel="Chemical shift (ppm)",
             ylabel="Intensity",
             xreversed=true,
             title="Overview")
    
    # Get the first slice
    x = expt.specdata.x[1]
    y = expt.specdata.y[1]
    
    # Plot spectrum
    lines!(ax, x, y; color=:black)
    
    # Add regions
    for (idx, region) in enumerate(expt.regions[])
        # Get region properties
        xstart = region.xstart[]
        xend = region.xend[]
        label = region.label[]
        color = region.color[]
        
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
            poly!(ax, Point2f.(xs, ys);
                 color=(color, 0.3),
                 strokewidth=2,
                 strokecolor=color)
            
            # Add label at top of region
            middle_x = (xstart + xend) / 2
            middle_y_idx = findfirst(xi -> xi >= middle_x, x_region)
            if !isnothing(middle_y_idx)
                middle_y = maximum(y_region) * 1.05
                text!(ax, middle_x, middle_y,
                     text=label,
                     color=color,
                     align=(:center, :bottom))
            end
        end
    end
    
    # Add baseline
    hlines!(ax, [0]; color=:gray, linestyle=:dash)
    
    # Save figure
    save(joinpath(folder, "overview.pdf"), fig)
    save(joinpath(folder, "overview.png"), fig; px_per_unit=2)
end

"""
    save_region_plot!(expt::Experiment1D, region::Region, folder::AbstractString)

Save a plot for a specific region.
"""
function save_region_plot!(expt::Experiment1D, region::Region, folder::AbstractString)
    fig = Figure(; size=(900, 600))
    
    # Create layout with spectrum and fit plots
    spectrum_panel = fig[1, 1] = GridLayout()
    fit_panel = fig[2, 1] = GridLayout()
    
    # Create spectrum plot
    ax_spectrum = Axis(spectrum_panel[1, 1];
                      xlabel="Chemical shift (ppm)",
                      ylabel="Intensity",
                      xreversed=true,
                      title="Region: $(region.label[])")
    
    # Get the first slice
    x = expt.specdata.x[1]
    y = expt.specdata.y[1]
    
    # Get region properties
    xstart = region.xstart[]
    xend = region.xend[]
    color = region.color[]
    
    # Find points within region
    idx_in_region = findall(xi -> xstart <= xi <= xend, x)
    
    if !isempty(idx_in_region)
        # Plot spectrum
        lines!(ax_spectrum, x, y; color=:black)
        
        # Create region fill
        x_region = x[idx_in_region]
        y_region = y[idx_in_region]
        
        # Create polygon vertices for fill
        xs = [x_region; reverse(x_region)]
        ys = [y_region; zeros(length(y_region))]
        
        # Add region visualization
        poly!(ax_spectrum, Point2f.(xs, ys);
             color=(color, 0.3),
             strokewidth=2,
             strokecolor=color)
        
        # Add baseline
        hlines!(ax_spectrum, [0]; color=:gray, linestyle=:dash)
        
        # Focus on region with some margin
        margin = (xend - xstart) * 0.2
        xlims!(ax_spectrum, (xend + margin, xstart - margin))
    end
    
    # Create experiment-specific fit plot
    create_fit_plot!(fig, fit_panel, expt, region)
    
    # Save figure
    safe_label = replace(region.label[], r"[^\w\s-]" => "_")
    save(joinpath(folder, "region_$(safe_label).pdf"), fig)
    save(joinpath(folder, "region_$(safe_label).png"), fig; px_per_unit=2)
end

"""
    create_fit_plot!(fig, panel, expt, region)

Create experiment-specific fit plot for a region.
Default implementation - can be overridden by specific experiment types.
"""
function create_fit_plot!(fig, panel, expt::Experiment1D, region::Region)
    # Create basic axis
    ax = Axis(panel[1, 1];
             xlabel="x",
             ylabel="y",
             title="Fit Results")
    
    # Add text with fit information
    if region.postfitted[]
        info_text = []
        
        # Add parameter info
        for (name, param) in region.postparameters
            value = param.value[][1]
            uncertainty = param.uncertainty[][1]
            push!(info_text, "$(param.label): $(round(value, digits=4)) ± $(round(uncertainty, digits=4))")
        end
        
        # Display text
        text!(ax, 0.05, 0.95;
             text=join(info_text, "\n"),
             align=(:left, :top),
             space=:relative)
    else
        text!(ax, 0.5, 0.5;
             text="No fit results available",
             align=(:center, :center),
             space=:relative)
    end
end

"""
    create_fit_plot!(fig, panel, expt::DiffusionExperiment, region::Region)

Create diffusion-specific fit plot.
"""
function create_fit_plot!(fig, panel, expt::DiffusionExperiment, region::Region)
    # Create axis
    ax = Axis(panel[1, 1];
             xlabel="Gradient strength / T m⁻¹",
             ylabel="Normalized integral",
             title="Diffusion Fit")
    
    if !region.postfitted[]
        text!(ax, 0.5, 0.5;
             text="No fit results available",
             align=(:center, :center),
             space=:relative)
        return
    end
    
    # Get data
    integrals = region.parameters[:integral].value[]
    errors = region.parameters[:error].value[]
    gradients = expt.gradients
    
    # Normalize by fitted amplitude
    amp = region.postparameters[:A].value[][1]
    integrals_norm = integrals ./ amp
    errors_norm = errors ./ amp
    
    # Create scatter plot of data
    scatter!(ax, gradients, integrals_norm;
            color=region.color[],
            markersize=8)
    
    # Add error bars
    for i in 1:length(gradients)
        lines!(ax, [gradients[i], gradients[i]], 
              [integrals_norm[i] - errors_norm[i], integrals_norm[i] + errors_norm[i]];
              color=region.color[])
    end
    
    # Create fit line
    g_range = range(0, maximum(gradients) * 1.1, 100)
    D = region.postparameters[:D].value[][1]
    y_fit = [exp(-(expt.γ * expt.δ * expt.σ * g)^2 * (expt.Δ - expt.δ/3) * D) for g in g_range]
    
    lines!(ax, g_range, y_fit;
          color=:red,
          linewidth=2)
    
    # Add fit results as text
    D_val = region.postparameters[:D].value[][1]
    D_err = region.postparameters[:D].uncertainty[][1]
    D_text = "D = $(round(D_val * 1e10, digits=2)) ± $(round(D_err * 1e10, digits=2)) × 10⁻¹⁰ m²/s"
    
    rH_text = ""
    if haskey(region.postparameters, :rH)
        rH_val = region.postparameters[:rH].value[][1]
        rH_err = region.postparameters[:rH].uncertainty[][1]
        rH_text = "\nrH = $(round(rH_val, digits=1)) ± $(round(rH_err, digits=1)) Å"
    end
    
    text!(ax, 0.05, 0.95;
         text=D_text * rH_text,
         align=(:left, :top),
         space=:relative)
end

"""
    create_fit_plot!(fig, panel, expt::TRACTExperiment, region::Region)

Create TRACT-specific fit plot.
"""
function create_fit_plot!(fig, panel, expt::TRACTExperiment, region::Region)
    # Create axis
    ax = Axis(panel[1, 1];
             xlabel="Relaxation time / s",
             ylabel="Normalized integral",
             title="TRACT Fit")
    
    if !region.postfitted[]
        text!(ax, 0.5, 0.5;
             text="No fit results available",
             align=(:center, :center),
             space=:relative)
        return
    end
    
    # Get data
    trosy_integrals = region.parameters[:integral_trosy].value[]
    trosy_errors = region.parameters[:error_trosy].value[]
    antitrosy_integrals = region.parameters[:integral_antitrosy].value[]
    antitrosy_errors = region.parameters[:error_antitrosy].value[]
    
    # Normalize by fitted amplitudes
    trosy_amp = region.postparameters[:A_trosy].value[][1]
    antitrosy_amp = region.postparameters[:A_antitrosy].value[][1]
    
    trosy_norm = trosy_integrals ./ trosy_amp
    trosy_err_norm = trosy_errors ./ trosy_amp
    antitrosy_norm = antitrosy_integrals ./ antitrosy_amp
    antitrosy_err_norm = antitrosy_errors ./ antitrosy_amp
    
    # Plot TROSY data
    scatter!(ax, expt.τ_trosy, trosy_norm;
            color=:blue,
            markersize=8,
            label="TROSY")
    
    # Plot anti-TROSY data
    scatter!(ax, expt.τ_antitrosy, antitrosy_norm;
            color=:red,
            markersize=8,
            label="Anti-TROSY")
    
    # Add error bars
    for i in 1:length(expt.τ_trosy)
        lines!(ax, [expt.τ_trosy[i], expt.τ_trosy[i]], 
              [trosy_norm[i] - trosy_err_norm[i], trosy_norm[i] + trosy_err_norm[i]];
              color=:blue)
    end
    
    for i in 1:length(expt.τ_antitrosy)
        lines!(ax, [expt.τ_antitrosy[i], expt.τ_antitrosy[i]], 
              [antitrosy_norm[i] - antitrosy_err_norm[i], antitrosy_norm[i] + antitrosy_err_norm[i]];
              color=:red)
    end
    
    # Create fit lines
    t_max = max(maximum(expt.τ_trosy), maximum(expt.τ_antitrosy))
    t_range = range(0, t_max * 1.1, 100)
    
    R2_trosy = region.postparameters[:R2_trosy].value[][1]
    R2_antitrosy = region.postparameters[:R2_antitrosy].value[][1]
    
    trosy_fit = [exp(-R2_trosy * t) for t in t_range]
    antitrosy_fit = [exp(-R2_antitrosy * t) for t in t_range]
    
    lines!(ax, t_range, trosy_fit;
          color=:blue,
          linewidth=2)
    
    lines!(ax, t_range, antitrosy_fit;
          color=:red,
          linewidth=2)
    
    # Add legend
    axislegend(ax; position=:best)
    
    # Add fit results as text
    R2_trosy_val = region.postparameters[:R2_trosy].value[][1]
    R2_trosy_err = region.postparameters[:R2_trosy].uncertainty[][1]
    R2_antitrosy_val = region.postparameters[:R2_antitrosy].value[][1]
    R2_antitrosy_err = region.postparameters[:R2_antitrosy].uncertainty[][1]
    τc_val = region.postparameters[:τc].value[][1]
    τc_err = region.postparameters[:τc].uncertainty[][1]
    
    result_text = """
    TROSY R₂: $(round(R2_trosy_val, digits=2)) ± $(round(R2_trosy_err, digits=2)) s⁻¹
    Anti-TROSY R₂: $(round(R2_antitrosy_val, digits=2)) ± $(round(R2_antitrosy_err, digits=2)) s⁻¹
    ΔR: $(round(region.postparameters[:ΔR].value[][1], digits=2)) ± $(round(region.postparameters[:ΔR].uncertainty[][1], digits=2)) s⁻¹
    τc: $(round(τc_val, digits=1)) ± $(round(τc_err, digits=1)) ns
    """
    
    text!(ax, 0.95, 0.95;
         text=result_text,
         align=(:right, :top),
         space=:relative)
end