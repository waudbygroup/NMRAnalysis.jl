# NMR Integration Functions

These functions provide tools for integrating regions of NMR spectra, with optional noise estimation capabilities.

## Main Function

```julia
integrate(filenames, shifts, width, noiseregion)
```
Integrates a region/regions of a 1D NMR spectrum or spectra centred at a specified chemical shift.
Computes an uncertainty if given a region of noise.
Integrals are scaled for receiver gain and number of scans.

### Arguments
- `filename`: Path to the NMR data file, or a list of filenames.
- `shift`: Chemical shift (in ppm) at the centre of the integration region, or a list of chemical shifts.
- `width`: Width of the integration region in ppm
- `noiseregion`: Optional tuple or range specifying a region (in ppm) for noise estimation

### Returns
- Without noise region: Returns the scaled integral value
- With noise region: Returns the integral with uncertainty (using the `±` operator)

### Details
The function performs the following operations:
1. Loads the NMR spectrum using `loadnmr`
2. Verifies the spectrum is 1-dimensional
3. Integrates a region of width `width` centred at `shift`
4. Scales the integral by the spectrum's scaling factor

When a noise region is provided, the function:
1. Divides the noise region into segments of the same width as the integration region
2. Calculates integrals for each noise segment
3. Computes the standard deviation of the noise integrals
4. Returns the main integral with uncertainty derived from the noise analysis

### Errors
- Throws an error if attempting to integrate multi-dimensional data

## Vector Methods

```julia
integrate(filename::String, shifts::Vector, width, noiseregion=nothing)
```

Integrates multiple regions in a single spectrum.

### Arguments
- `shifts::Vector`: Vector of chemical shifts (in ppm) for multiple integration regions
- Other arguments as per the main function

### Returns
- Vector of integral values (with or without uncertainties)

```julia
integrate(filenames::Vector{String}, shifts, width, noiseregion=nothing)
```

Integrates regions across multiple spectra.

### Arguments
- `filenames::Vector{String}`: Vector of paths to NMR data files
- Other arguments as per the main function

### Returns
- Vector of integral values (with or without uncertainties)

## Examples

```julia
# Simple integration
I = integrate("1d_spectrum.nmr", 7.24, 0.1)

# Integration with noise estimation
I = integrate("1d_spectrum.nmr", 7.24, 0.1, (0.0..1.0))

# Multiple regions in one spectrum
shifts = [7.24, 3.72, 2.1]
Is = integrate("1d_spectrum.nmr", shifts, 0.1, (0.0..1.0))

# Same region across multiple spectra
spectra = ["spectrum1.nmr", "spectrum2.nmr", "spectrum3.nmr"]
Is = integrate(spectra, 7.24, 0.1, (0.0..1.0))
```

## Notes
- The noise estimation uses equally spaced integration regions of the same width as the signal region
- Integration results are automatically scaled using the spectrum's scaling factor
- When using noise estimation, results are returned with uncertainty using the `±` operator
- The function assumes the input spectra are properly phased and baseline corrected