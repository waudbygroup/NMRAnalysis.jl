"""
    lpre2d(inputfilenames, paramagnetic_concs, expttype, Trelax)

Create an LPRE experiment from files and experimental parameters and launch analysis interface.
`expttype` should be `:hsqc` or `:hmqc`. `Trelax` is the timing during which relaxation
can occur during the sequence (magnetisation transfer delays etc) - this is specific to
the pulse sequence used.

Can be used to analyse solvent PREs, in which case concentrations should be specified -
or to analyse protein PREs, in which case concentrations should be set to 0 and 1 for
diamagnetic and paramagnetic states respectively.
"""
function lpre2d(inputfilenames, paramagnetic_concs, recycle_delays, expttype, Trelax)
    expt = LPREExperiment(inputfilenames, paramagnetic_concs, recycle_delays, expttype, Trelax)
    gui!(expt)
end

"""
    LPREExperiment

Longitudinal paramagnetic relaxation enhancement experiment.

# Fields
- `specdata`: Spectral data and metadata
- `peaks`: Observable list of peaks  
- `paramagnetic_concs`: Vector of paramagnetic agent concentrations
- `recycle_delays`: Vector of recycle delays
- `expttype`: Either :hsqc or :hmqc
- `Trelax`: Relaxation time
"""
struct LPREExperiment <: FixedPeakExperiment
    specdata::Any
    peaks::Any
    paramagnetic_concs::Any
    recycle_delays::Any
    expttype::Any
    Trelax::Any

    clusters::Any
    touched::Any
    isfitting::Any

    xradius::Any
    yradius::Any
    state::Any

    function LPREExperiment(specdata, peaks, paramagnetic_concs, recycle_delays, expttype, Trelax)
        expt = new(specdata, peaks, paramagnetic_concs, recycle_delays, expttype, Trelax,
                   Observable(Vector{Vector{Int}}()), # clusters
                   Observable(Vector{Bool}()), # touched
                   Observable(true), # isfitting
                   Observable(0.03; ignore_equal_values=true), # xradius
                   Observable(0.2; ignore_equal_values=true), # yradius
                   Observable{Dict}())
        setupexptobservables!(expt)
        expt.state[] = preparestate(expt)
        expt
    end
end

"""
    LPREExperiment(inputfilenames, paramagnetic_concs, recycle_delays, expttype, Trelax)

Create a LPRE experiment from files and experimental parameters.
`expttype` should be `:hsqc` or `:hmqc`. `Trelax` is the timing during which relaxation
can occur during the sequence (magnetisation transfer delays etc) - this is specific to
the pulse sequence used.

Can be used to analyse solvent PREs, in which case concentrations should be specified -
or to analyse protein PREs, in which case concentrations should be set to 0 and 1 for
diamagnetic and paramagnetic states respectively.
"""
function LPREExperiment(inputfilenames, paramagnetic_concs, recycle_delays, exptexperimenttype, Trelax)
    exptexperimenttype in [:hsqc, :hmqc] ||
        throw(ArgumentError("Experiment type must be :hsqc or :hmqc"))

    specdata = preparespecdata(inputfilenames, paramagnetic_concs, recycle_delays, LPREExperiment)
    peaks = Observable(Vector{Peak}())

    return LPREExperiment(specdata, peaks, paramagnetic_concs, recycle_delays, exptexperimenttype, Trelax)
end

# load the NMR data and prepare the SpecData object
function preparespecdata(inputfilenames, paramagnetic_concs, recycle_delays, ::Type{LPREExperiment})
    @debug "Preparing spec data for LPRE relaxation experiment: $inputfilenames"

    spectra = map(loadnmr, inputfilenames)
    x = map(spec -> data(spec, F1Dim), spectra)
    y = map(spec -> data(spec, F2Dim), spectra)
    z = map(spec -> data(spec) / scale(spec), spectra)
    σ = map(spec -> spec[:noise] / scale(spec), spectra)

    zlabels = map(zip(paramagnetic_concs, recycle_delays)) do (conc, d)
        if conc == 0
            "Diamagnetic (d1 = $(d)s)"
        else 
            "Paramagnetic (conc = $conc, d1 = $(d)s)"
        end
    end

    return SpecData(spectra, x, y,
                    z ./ σ[1],
                    σ ./ σ[1],
                    zlabels)
end

"""Add peak to experiment, setting up type-specific parameters."""
function addpeak!(expt::LPREExperiment, initialposition::Point2f, label="",
                  xradius=expt.xradius[], yradius=expt.yradius[])
    expt.state[][:total_peaks][] += 1
    if label == ""
        label = "X$(expt.state[][:total_peaks][])"
    end
    @debug "Add peak $label at $initialposition"
    newpeak = Peak(initialposition, label, xradius, yradius)

    # pars: R1x, R2x, R2y, amp
    R1x0 = MaybeVector(2.0)
    R2x0 = MaybeVector(30.0)
    R2y0 = MaybeVector(15.0)
    R1x = Parameter("R1x", R1x0; minvalue=0.2, maxvalue=100.0)
    R2x = Parameter("R2x", R2x0; minvalue=1.0, maxvalue=100.0)
    R2y = Parameter("R2y", R2y0; minvalue=1.0, maxvalue=100.0)

    # get initial values for amplitude
    x0, y0 = initialposition
    ix = findnearest(expt.specdata.x[1], x0)
    iy = findnearest(expt.specdata.y[1], y0)
    amp0 = expt.specdata.z[1][ix, iy]
    amp0 = MaybeVector(amp0)
    amp = Parameter("Amplitude", amp0)

    # PRE
    Γ1 = Parameter("Γ1", MaybeVector(2.0); minvalue=0.0, maxvalue=200.0)
    Γ2 = Parameter("Γ2", MaybeVector(10.0); minvalue=0.0, maxvalue=200.0)

    newpeak.parameters[:R1x] = R1x
    newpeak.parameters[:R2x] = R2x
    newpeak.parameters[:R2y] = R2y
    newpeak.parameters[:amp] = amp
    newpeak.parameters[:Γ1] = Γ1
    newpeak.parameters[:Γ2] = Γ2

    newpeak.postparameters[:Γ1] = Parameter("Γ1", 0.0)
    newpeak.postparameters[:Γ2] = Parameter("Γ2", 0.0)

    push!(expt.peaks[], newpeak)
    return notify(expt.peaks)
end

"""Simulate single peak according to experiment type."""
function simulate!(z, peak::Peak, expt::LPREExperiment, xbounds=nothing, ybounds=nothing)
    R1x0 = peak.parameters[:R1x].value[][1]
    R2x0 = peak.parameters[:R2x].value[][1]
    R2y0 = peak.parameters[:R2y].value[][1]
    amp0 = peak.parameters[:amp].value[][1]
    Γ1 = peak.parameters[:Γ1].value[][1]
    Γ2 = peak.parameters[:Γ2].value[][1]

    for i in 1:nslices(expt)
        # get axis references for window functions
        xaxis = dims(expt.specdata.nmrdata[i], F1Dim)
        yaxis = dims(expt.specdata.nmrdata[i], F2Dim)
        # get axis shift values
        x = isnothing(xbounds) ? expt.specdata.x[i] : expt.specdata.x[i][xbounds[i]]
        y = isnothing(ybounds) ? expt.specdata.y[i] : expt.specdata.y[i][ybounds[i]]

        x0 = peak.parameters[:x].value[][i]
        y0 = peak.parameters[:y].value[][i]

        # apply PRE to linewidths and amplitude
        R1x = R1x0 + Γ1 * expt.paramagnetic_concs[i]
        R2x = R2x0 + Γ2 * expt.paramagnetic_concs[i]
        if expt.expttype == :hmqc
            R2y = R2y0 + Γ2 * expt.paramagnetic_concs[i]
        else
            R2y = R2y0 #+ Γ1 * expt.paramagnetic_concs[i] # TODO - check this?
        end
        amp = amp0 * exp(-Γ2 * expt.paramagnetic_concs[i] * expt.Trelax)
        amp *= (1 - exp(-R1x * expt.recycle_delays[i]))

        # find indices of x and y axes within peak radius of peak position
        xi = x0 .- peak.xradius[] .≤ x .≤ x0 .+ peak.xradius[]
        yi = y0 .- peak.yradius[] .≤ y .≤ y0 .+ peak.yradius[]
        xs = x[xi]
        ys = y[yi]
        # NB. scale intensities by R2x and R2y to decouple amplitude estimation from linewidth
        zx = NMRTools.NMRBase._lineshape(getω(xaxis, x0), R2x, getω(xaxis, xs),
                                         xaxis[:window], RealLineshape())
        zy = (π^2 * amp * R2x0 * R2y0) *
             NMRTools.NMRBase._lineshape(getω(yaxis, y0), R2y, getω(yaxis, ys),
                                         yaxis[:window], RealLineshape())
        z[i][xi, yi] .+= zx .* zy'
    end
end

"""Calculate final parameters after fitting."""
function postfit!(peak::Peak, ::LPREExperiment)
    peak.postparameters[:Γ1].uncertainty[] .= peak.parameters[:Γ1].uncertainty[]
    peak.postparameters[:Γ1].value[] .= peak.parameters[:Γ1].value[]
    peak.postparameters[:Γ2].uncertainty[] .= peak.parameters[:Γ2].uncertainty[]
    peak.postparameters[:Γ2].value[] .= peak.parameters[:Γ2].value[]
    peak.postfitted[] = true
end

"""Return descriptive text for slice idx."""
function slicelabel(expt::LPREExperiment, idx)
    return "$(expt.specdata.zlabels[idx]) ($idx of $(nslices(expt)))"
end

"""Return formatted text describing peak idx."""
function peakinfotext(expt::LPREExperiment, idx)
    if idx == 0
        return "No peak selected"
    end
    peak = expt.peaks[][idx]
    if peak.postfitted[]
        return "Peak: $(peak.label[])\n" *
               "Γ1: $(peak.parameters[:Γ1].value[][1] ± peak.parameters[:Γ1].uncertainty[][1]) s⁻¹ [conc⁻¹]\n" *
               "Γ2: $(peak.parameters[:Γ2].value[][1] ± peak.parameters[:Γ2].uncertainty[][1]) s⁻¹ [conc⁻¹]\n" *
               "\n" *
               "δX: $(peak.parameters[:x].value[][1] ± peak.parameters[:x].uncertainty[][1]) ppm\n" *
               "δY: $(peak.parameters[:y].value[][1] ± peak.parameters[:y].uncertainty[][1]) ppm\n" *
               "Amplitude: $(peak.parameters[:amp].value[][1] ± peak.parameters[:amp].uncertainty[][1])\n" *
               "X Linewidth: $(peak.parameters[:R2x].value[][1] ± peak.parameters[:R2x].uncertainty[][1]) s⁻¹\n" *
               "Y Linewidth: $(peak.parameters[:R2y].value[][1] ± peak.parameters[:R2y].uncertainty[][1]) s⁻¹"
    else
        return "Peak: $(peak.label[])\n" *
               "Not fitted"
    end
end

"""Return formatted text describing experiment."""
function experimentinfo(expt::LPREExperiment)
    return "Analysis type: PRE experiment\n" *
           "Filename: $(expt.specdata.nmrdata[1][:filename])\n" *
           "PRE agent concs: $(join(expt.paramagnetic_concs, ", "))\n" *
           "Recycle delays: $(join(expt.recycle_delays, ", ")) s\n" *
           "Experiment type: $(expt.expttype==:hsqc ? "HSQC" : "HMQC")\n" *
           "Relaxation time: $(expt.Trelax)\n" *
           "Number of peaks: $(length(expt.peaks[]))\n" *
           "Experiment title: $(expt.specdata.nmrdata[1][:title])\n"
end


##
cd("/Users/chris/NMR/crick-800/chris_hewl_240429/")
recovery2d(["103","104","105","106","107","108","109","110","111"],[.1,.2,.3,.4,.5,0.75,1,1.25,1.5])
recovery2d(["112","113","114","115","116","117","118","119","120"],[.1,.2,.3,.4,.5,0.75,1,1.25,1.5])
recovery2d(["10","12","14","16","18","20","22","24","26"],[.1,.2,.3,.4,.5,0.75,1,1.25,1.5])
recovery2d(["11","13","15","17","19","21","23","25","27"],[.1,.2,.3,.4,.5,0.75,1,1.25,1.5])
lpre2d(["10","12","14","16","18","20","22","24","26", "103","104","105","106","107","108","109","110","111"],
           [2, 2, 2, 2, 2, 2, 2, 2, 2, 0.,0,0,0,0,0,0,0,0],
           [0.1,0.2,0.3,0.4,0.5,0.75,1,1.25,1.5, 0.1,0.2,0.3,0.4,0.5,0.75,1,1.25,1.5],
           :hmqc, 0.01)
lpre2d(["11","13","15","17","19","21","23","25","27", "112","113","114","115","116","117","118","119","120"],
           [2, 2, 2, 2, 2, 2, 2, 2, 2, 0.,0,0,0,0,0,0,0,0],
           [0.1,0.2,0.3,0.4,0.5,0.75,1,1.25,1.5, 0.1,0.2,0.3,0.4,0.5,0.75,1,1.25,1.5],
           :hmqc, 0.01)
