"""
    integrate(filename::String, shift::Number, width, noiseregion=nothing)
    integrate(filename::String, shifts::Vector, width, noiseregion=nothing)
    integrate(filenames::Vector{String}, shifts, width, noiseregion=nothing)

Integrate regions of 1D NMR spectra with optional noise estimation.

# Arguments
- `filename::String` or `filenames::Vector{String}`: Path(s) to NMR data file(s)
- `shift::Number` or `shifts::Vector`: Chemical shift(s) in ppm at the centre of integration region(s)
- `width`: Width of integration region in ppm
- `noiseregion=nothing`: Optional range for noise estimation (e.g., `0.0..1.0`)

# Returns
Without noise region:
- Single shift: Scaled integral value
- Multiple shifts/files: Vector of scaled integral values

With noise region:
- Single shift: Integral ± uncertainty
- Multiple shifts/files: Vector of integral ± uncertainty values

# Examples
```julia
# Single region integration
I = integrate("spectrum.nmr", 7.24, 0.1)

# With noise estimation
I = integrate("spectrum.nmr", 7.24, 0.1, (0.0..1.0))

# Multiple regions
Is = integrate("spectrum.nmr", [7.24, 3.72], 0.1, (0.0..1.0))

# Multiple spectra
Is = integrate(["spec1.nmr", "spec2.nmr"], 7.24, 0.1)
```

# Notes
- Only works with 1D spectra
- Integration regions are centred at specified chemical shifts
- Noise estimation uses equally spaced regions of width `width`
- Results are automatically scaled by the spectrum's scaling factor
- Assumes properly phased and baseline-corrected spectra

See also: [`loadnmr`](@ref)
"""
function integrate(filename::String, shift::Number, width::Number, noiseregion::Union{Nothing,Tuple}=nothing)
    spec = loadnmr(filename)
    ndims(spec) == 1 || error("attempting to integrate $(ndims(spec))D experiment: $filename")

    reg = shift - 0.5width .. shift + 0.5width
    I = sum(spec[reg]) / scale(spec)

    if !isnothing(noiseregion)
        # get points equally spaced by width within the noise region
        n1, n2 = minmax(noiseregion...)
        noisepoints = (n1+0.5width):width:n2
        noiseregions = map(x -> x-0.5width .. x+0.5width, noisepoints)
        Inoise = map(x->sum(spec[x]), noiseregions) ./ scale(spec)
        σ = std(Inoise)
        I = I ± σ
    end
    I /= 1e6

    println("$(basename(filename)):\t $(round(shift, digits=3)) ppm, $I")

    return (shift, I)
end

function integrate(filename::String, width::Number, noiseregion::Union{Nothing,Tuple}=nothing; SNR_threshold=8)
    spec = loadnmr(filename)
    px = findpeaks1d(spec, SNR_threshold)
    integrate(filename, px, width, noiseregion)
end

function integrate(filename::String, shifts::Vector, width::Number, noiseregion::Union{Nothing,Tuple}=nothing)
    map(s -> integrate(filename, s, width, noiseregion), shifts)
end

function integrate(filenames::Vector{String}, shifts::Vector, width::Number, noiseregion::Union{Nothing,Tuple}=nothing)
    map(f -> integrate(f, shifts, width, noiseregion), filenames)
end

function integrate(filenames::Vector{String}, width::Number, noiseregion::Union{Nothing,Tuple}=nothing; SNR_threshold=8)
    map(f -> integrate(f, width, noiseregion; SNR_threshold=SNR_threshold), filenames)
end

function findpeaks1d(spec, snr_threshold=5)
    x = data(spec, F1Dim)
    y = data(spec[:,1]) / spec[:noise]

    pks = findmaxima(y)     # detect peaks
    peakheights!(pks; min=snr_threshold) # filter list

    return x[pks.indices]
end