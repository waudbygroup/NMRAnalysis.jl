struct RelaxationExperiment <: Experiment
    specdata
    peaks
    relaxationtimes

    adjacancy
    clusters
    RelaxationExperiment(specdata, peaks, relaxationtimes) = setup_graph_fields(new(specdata, peaks, relaxationtimes))
end

export RelaxationExperiment

# create a new relaxation experiment
function RelaxationExperiment(inputfilename, taufilename)
    relaxationtimes = vec(readdlm(taufilename, comments=true))
    specdata = preparespecdata(inputfilename, relaxationtimes, ::RelaxationExperiment)
    peaks = Observable(Vector{Peak}())

    RelaxationExperiment(specdata, peaks, relaxationtimes)
end


# implementation requirements
hasfixedpositions(expt::RelaxationExperiment) = true

function addpeak!(expt::RelaxationExperiment, position)
    # TODO set up initial parameters
    # push!(expt.peaks[], Peak(position))
end

function simulate!(z, x, y, peak, expt::RelaxationExperiment)
    n = length(z)
    xaxis = data(expt.specdata.nmrdata[1], F1Dim)
    yaxis = data(expt.specdata.nmrdata[1], F2Dim)
    
    for i in 1:n
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

function mask!(z, x, y, peak, expt::RelaxationExperiment)
    n = length(z)
    for i in 1:n
        maskellipse!(z[i], x[i], y[i],
            peak.initialposition[][1],
            peak.initialposition[][2],
            peak.xradius[], peak.yradius[])
    end
end

# load the NMR data and prepare the SpecData object
function preparespecdata(inputfilename, relaxationtimes, ::Type{RelaxationExperiment})
    spec = loadnmr(inputfilename)
    x = SingleElementVector(data(spec, F1Dim))
    y = SingleElementVector(data(spec, F2Dim))
    
    dat = data(spec) / scale(spec)
    σ = spec[:noise] / scale(spec)

    z = eachslice(dat, dims=3)
    zlabels = map(t -> "τ = $t", relaxationtimes)

    SpecData(SingleElementVector(spec), x, y, z, σ, zlabels)
end