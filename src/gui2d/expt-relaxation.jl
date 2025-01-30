struct RelaxationExperiment <: Experiment
    specdata
    peaks
    relaxationtimes

    clusters
    touched
    isfitting

    RelaxationExperiment(specdata, peaks, relaxationtimes) = begin
        expt = new(specdata, peaks, relaxationtimes,
            Observable(Vector{Vector{Int}}()),
            Observable(Vector{Bool}()),
            Observable(false)
            )
        setupexptobservables!(expt)
        expt
    end
end

export RelaxationExperiment

# create a new relaxation experiment
function RelaxationExperiment(inputfilename, taufilename)
    relaxationtimes = vec(readdlm(taufilename, comments=true))
    specdata = preparespecdata(inputfilename, relaxationtimes, RelaxationExperiment)
    peaks = Observable(Vector{Peak}())

    RelaxationExperiment(specdata, peaks, relaxationtimes)
end


# implementation requirements
hasfixedpositions(expt::RelaxationExperiment) = true

function addpeak!(expt::RelaxationExperiment, initialposition::Point2f, label, xradius=0.03, yradius=0.3)
    newpeak = Peak(initialposition, label, xradius, yradius)
    # pars: R2x, R2y, amp
    R2x0 = MaybeVector(10.)
    R2y0 = MaybeVector(10.)
    R2x = Parameter("R2x", R2x0, 1., 100.)
    R2y = Parameter("R2y", R2y0, 1., 100.)
    # get initial values for amplitude
    x0, y0 = initialposition
    amp0 = map(1:nslices(expt)) do i
        ix = findnearest(expt.specdata.x[i], x0)
        iy = findnearest(expt.specdata.y[i], y0)
        expt.specdata.z[i][ix, iy] / (R2x0[i] * R2y0[i])
    end
    amp = Parameter("amp", amp0, -1e4, 1e4)
    newpeak.parameters[:R2x] = R2x
    newpeak.parameters[:R2y] = R2y
    newpeak.parameters[:amp] = amp

    @debug "Adding peak $label at $initialposition" newpeak
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
        zx = NMRTools.NMRBase._lineshape(getω(xaxis, x0), R2x, getω(xaxis, xs), xaxis[:window], RealLineshape())
        zy = amp * NMRTools.NMRBase._lineshape(getω(yaxis, y0), R2y, getω(yaxis, ys), yaxis[:window], RealLineshape())
        z[i][xi, yi] .+= zx .* zy'
    end
end

function mask!(z, peak::Peak, expt::RelaxationExperiment)
    @debug "masking peak $peak"
    n = length(z)
    for i in 1:n
        @show i
        x = data(expt.specdata.nmrdata[i], F1Dim)
        y = data(expt.specdata.nmrdata[i], F2Dim)
        maskellipse!(z[i], x, y,
            initialposition(peak)[][i][1],
            initialposition(peak)[][i][2],
            peak.xradius[], peak.yradius[])
    end
end

# load the NMR data and prepare the SpecData object
function preparespecdata(inputfilename, relaxationtimes, ::Type{RelaxationExperiment})
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