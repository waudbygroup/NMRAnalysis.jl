
function setup_expt_fields(expt)
    expt.adjacency = @lift makeadjacency($(expt.peaks), expt)
    expt.clusters = @lift findclusters($(expt.adjacency))

    expt
end

mutable struct RelaxationExperiment <: Experiment
    specdata
    peaks
    relaxationtimes

    adjacency
    clusters
    RelaxationExperiment(specdata, peaks, relaxationtimes) = setup_expt_fields(new(specdata, peaks, relaxationtimes))
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
    R2x0 = 10.
    R2y0 = 10.
    R2x = Parameter("R2x", R2x0, 1., 100.)
    R2y = Parameter("R2y", R2y0, 1., 100.)
    # get initial values for amplitude
    x0, y0 = initialposition
    amp0 = map(1:nslices(expt)) do i
        ix = findnearest(expt.specdata.x[i], x0)
        iy = findnearest(expt.specdata.y[i], y0)
        expt.specdata.z[i][ix, iy] / (R2x0 * R2y0)
    end
    amp = Parameter("amp", amp0, -1e4, 1e4)
    newpeak.parameters[:R2x] = R2x
    newpeak.parameters[:R2y] = R2y
    newpeak.parameters[:amp] = amp

    @debug "Adding peak $label at $initialposition" newpeak
    push!(expt.peaks[], newpeak)
end

function simulate!(z, peak, expt::RelaxationExperiment)
    n = length(z)
    for i in 1:n
        # get axis references for window functions
        xaxis = dims(expt.specdata.nmrdata[i], F1Dim)
        yaxis = dims(expt.specdata.nmrdata[i], F2Dim)
        # get axis shift values
        x = data(xaxis)
        y = data(yaxis)

        pos = peak.initialposition[][i]
        x0 = pos[1]
        y0 = pos[2]
        R2x = peak.pars[:R2x][i].value
        R2y = peak.pars[:R2y][i].value
        amp = peak.pars[:amp][i].value
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
    @info "masking peak $peak"
    n = length(z)
    for i in 1:n
        x = data(expt.specdata.nmrdata[i], F1Dim)
        y = data(expt.specdata.nmrdata[i], F2Dim)
        maskellipse!(z[i], x, y,
            peak.initialposition[][i][1],
            peak.initialposition[][i][2],
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