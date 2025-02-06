struct RelaxationExperiment <: FixedPeakExperiment
    specdata
    peaks
    relaxationtimes

    clusters
    touched
    isfitting

    xradius
    yradius
    state

    RelaxationExperiment(specdata, peaks, relaxationtimes) = begin
        expt = new(specdata, peaks, relaxationtimes,
            Observable(Vector{Vector{Int}}()), # clusters
            Observable(Vector{Bool}()), # touched
            Observable(true), # isfitting
            Observable(0.03, ignore_equal_values=true), # xradius
            Observable(0.2, ignore_equal_values=true), # yradius
            Observable{Dict}()
            )
        setupexptobservables!(expt)
        expt.state[] = preparestate(expt)
        expt
    end
end

# create a new relaxation experiment
function RelaxationExperiment(inputfilename, taufilename)
    relaxationtimes = vec(readdlm(taufilename, comments=true))
    specdata = preparespecdata(inputfilename, relaxationtimes, RelaxationExperiment)
    peaks = Observable(Vector{Peak}())

    RelaxationExperiment(specdata, peaks, relaxationtimes)
end



# load the NMR data and prepare the SpecData object
function preparespecdata(inputfilename, relaxationtimes, ::Type{RelaxationExperiment})
    @debug "Preparing spec data for relaxation experiment: $inputfilename"
    spec = loadnmr(inputfilename)
    x = data(spec, F1Dim)
    y = data(spec, F2Dim)
    
    dat = data(spec) / scale(spec)
    σ = spec[:noise] / scale(spec)

    z = eachslice(dat, dims=3)
    zlabels = map(t -> "τ = $t", relaxationtimes)

    SpecData(SingleElementVector(spec),
        SingleElementVector(x),
        SingleElementVector(y),
        z ./ σ,
        SingleElementVector(1),
        zlabels)
end




# implementation requirements
function addpeak!(expt::RelaxationExperiment, initialposition::Point2f, label="",
                    xradius=expt.xradius[], yradius=expt.yradius[])
    expt.state[][:total_peaks][] += 1
    if label == ""
        label = "X$(expt.state[][:total_peaks][])"
    end
    @debug "Add peak $label at $initialposition"
    newpeak = Peak(initialposition, label, xradius, yradius)
    # pars: R2x, R2y, amp
    R2x0 = MaybeVector(10.)
    R2y0 = MaybeVector(10.)
    R2x = Parameter("R2x", R2x0, minvalue=1., maxvalue=100.)
    R2y = Parameter("R2y", R2y0, minvalue=1., maxvalue=100.)
    # get initial values for amplitude
    x0, y0 = initialposition
    amp0 = map(1:nslices(expt)) do i
        ix = findnearest(expt.specdata.x[i], x0)
        iy = findnearest(expt.specdata.y[i], y0)
        expt.specdata.z[i][ix, iy]
    end
    amp = Parameter("Amplitude", amp0)
    newpeak.parameters[:R2x] = R2x
    newpeak.parameters[:R2y] = R2y
    newpeak.parameters[:amp] = amp

    newpeak.postparameters[:relaxationrate] = Parameter("Relaxation rate", 4/maximum(expt.relaxationtimes))
    newpeak.postparameters[:amp] = Parameter("Amplitude", maximum(amp0))

    push!(expt.peaks[], newpeak)
    notify(expt.peaks)
end

function simulate!(z, peak::Peak, expt::RelaxationExperiment)
    n = length(z)
    for i in 1:n
        # get axis references for window functions
        xaxis = dims(expt.specdata.nmrdata[i], F1Dim)
        yaxis = dims(expt.specdata.nmrdata[i], F2Dim)
        # get axis shift values
        x = data(xaxis)
        y = data(yaxis)

        x0 = peak.parameters[:x].value[][i]
        y0 = peak.parameters[:y].value[][i]
        R2x = peak.parameters[:R2x].value[][i]
        R2y = peak.parameters[:R2y].value[][i]
        amp = peak.parameters[:amp].value[][i]
        # find indices of x and y axes within peak radius of peak position
        xi = x0 .- peak.xradius[] .≤ x .≤ x0 .+ peak.xradius[]
        yi = y0 .- peak.yradius[] .≤ y .≤ y0 .+ peak.yradius[]
        xs = x[xi]
        ys = y[yi]
        # NB. scale intensities by R2x and R2y to decouple amplitude estimation from linewidth
        zx = NMRTools.NMRBase._lineshape(getω(xaxis, x0), R2x, getω(xaxis, xs), xaxis[:window], RealLineshape())
        zy = (amp * R2x * R2y) * NMRTools.NMRBase._lineshape(getω(yaxis, y0), R2y, getω(yaxis, ys), yaxis[:window], RealLineshape())
        z[i][xi, yi] .+= zx .* zy'
    end
end



function postfit!(peak::Peak, expt::RelaxationExperiment)
    @debug "Post-fitting peak $(peak.label)" maxlog=10
    t = expt.relaxationtimes
    y = peak.parameters[:amp].value[]

    A0 = maximum(y)
    peak.postparameters[:amp].initialvalue[] .= A0
    R0 = peak.postparameters[:relaxationrate].initialvalue[][1]

    model(t, p) = p[1] * exp.(-p[2] * t)
    p0 = [A0, R0]
    fit = curve_fit(model, t, y, p0)
    pfit = coef(fit)
    perr = stderror(fit)

    peak.postparameters[:amp].value[] .= pfit[1]
    peak.postparameters[:amp].uncertainty[] .= perr[1]
    peak.postparameters[:relaxationrate].value[] .= pfit[2]
    peak.postparameters[:relaxationrate].uncertainty[] .= perr[2]

    peak.postfitted[] = true
end

function slicelabel(expt::RelaxationExperiment, idx)
    "τ = $(round(expt.relaxationtimes[idx],sigdigits=3)) s ($idx of $(nslices(expt)))"
end

function peakinfotext(expt::RelaxationExperiment, idx)
    if idx == 0
        return "No peak selected"
    end
    peak = expt.peaks[][idx]
    if peak.postfitted[]
        return "Peak: $(peak.label[])\n" *
            "Relaxation rate: $(peak.postparameters[:relaxationrate].value[][1] ± peak.postparameters[:relaxationrate].uncertainty[][1]) s⁻¹\n" *
            "\n" *
            "δX: $(peak.parameters[:x].value[][1] ± peak.parameters[:x].uncertainty[][1]) ppm\n" *
            "δY: $(peak.parameters[:y].value[][1] ± peak.parameters[:y].uncertainty[][1]) ppm\n" *
            "Amplitude: $(peak.postparameters[:amp].value[][1] ± peak.postparameters[:amp].uncertainty[][1])\n" *
            "X Linewidth: $(peak.parameters[:R2x].value[][1] ± peak.parameters[:R2x].uncertainty[][1]) s⁻¹\n" *
            "Y Linewidth: $(peak.parameters[:R2y].value[][1] ± peak.parameters[:R2y].uncertainty[][1]) s⁻¹"
    else
        return "Peak: $(peak.label[])\n" *
            "Not fitted"
    end
end



function experimentinfo(expt::RelaxationExperiment)
    "Analysis type: Relaxation experiment\n" *
    "Filename: $(expt.specdata.nmrdata[1][:filename])\n" *
    "Relaxation times: $(join(expt.relaxationtimes, ", ")) s\n" *
    "Number of peaks: $(length(expt.peaks[]))\n" *
    "Experiment title: $(expt.specdata.nmrdata[1][:title])\n"
end



function completestate!(state, expt::RelaxationExperiment)
    # Set up observables for the GUI
    state[:peakplot_obs_xy] = lift(state[:current_peak]) do peak
        isnothing(peak) ? Point2f[] : peak_plot_data(peak, expt)[1]
    end
    
    state[:peakplot_obs_xye] = lift(state[:current_peak]) do peak
        isnothing(peak) ? [(0.,0.,0.)] : peak_plot_data(peak, expt)[2]
    end
    
    state[:peakplot_fit_xy] = lift(state[:current_peak]) do peak
        isnothing(peak) ? Point2f[] : peak_plot_data(peak, expt)[3]
    end
end


"""
    peak_plot_data(peak, expt::RelaxationExperiment)

Extract plotting data for a single peak, returning observed points, error bars,
and fit line data.
"""
function peak_plot_data(peak, expt::RelaxationExperiment)
    # Calculate observed data points with errors
    t = expt.relaxationtimes
    y = peak.parameters[:amp].value[]
    err = peak.parameters[:amp].uncertainty[]
    obs_points = Point2f.(t, y)
    obs_errors = [(t[i], y[i], err[i]) for i in 1:length(t)]
    
    # Calculate fit line
    tpred = range(0, 1.1 * maximum(t), 100)
    A = peak.postparameters[:amp].value[][1]
    R = peak.postparameters[:relaxationrate].value[][1]
    ypred = A * exp.(-R * tpred)
    fit_points = Point2f.(tpred, ypred)
    
    return (obs_points, obs_errors, fit_points)
end


"""
    plot_peak!(ax, peak, relaxationtimes)

Plot a single peak's data and fit onto the given axis.
"""
function plot_peak!(ax, peak, expt::RelaxationExperiment)
    obs_points, obs_errors, fit_points = peak_plot_data(peak, expt)
    
    hlines!(ax, [0], linewidth=0)
    lines!(ax, fit_points, label="Fit", color=:red)
    errorbars!(ax, obs_errors, whiskerwidth=10)
    scatter!(ax, obs_points, label="Observed")
end


"""
    makepeakplot!(gui, state, expt::RelaxationExperiment)

Create interactive peak plot in GUI context.
"""
function makepeakplot!(gui, state, expt::RelaxationExperiment)
    gui[:axpeakplot] = ax = Axis(gui[:panelpeakplot][1,1], 
                                xlabel="Relaxation time / s", 
                                ylabel="Amplitude")
    
    hlines!(ax, [0], linewidth=0)
    lines!(ax, state[:peakplot_fit_xy], label="Fit", color=:red)
    errorbars!(ax, state[:peakplot_obs_xye], whiskerwidth=10)
    scatter!(ax, state[:peakplot_obs_xy], label="Observed")
end


"""
    save_peak_plots(expt::RelaxationExperiment, folder::AbstractString)

Save individual plots for each peak in the experiment to the specified folder.
"""
function save_peak_plots!(expt::RelaxationExperiment, folder::AbstractString)
    CairoMakie.activate!()
    for peak in expt.peaks[]
        fig = Figure()
        ax = Axis(fig[1,1], 
                  xlabel="Relaxation time / s", 
                  ylabel="Amplitude",
                  title="$(peak.label[])")
        
        plot_peak!(ax, peak, expt)
        # axislegend(ax)
        
        save(joinpath(folder, "peak_$(peak.label[]).pdf"), fig)
    end
    GLMakie.activate!()
end