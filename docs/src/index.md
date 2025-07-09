# NMRAnalysis.jl

A Julia package for analysis of NMR data, with particular focus on relaxation, diffusion, screening and
protein dynamics experiments.

!!! note
    This package is under active development and its features and syntax may change.

## Overview

NMRAnalysis.jl provides tools for processing and analyzing various types of NMR experiments commonly
used in biomolecular NMR:

- Relaxation measurements (T1, T2)
- Heteronuclear NOE
- Diffusion
- TRACT (TROSY/anti-TROSY Cross-Correlated Relaxation)

## Quick Start

```julia
using NMRAnalysis

# Set your working directory to where your NMR data is located
cd("/path/to/nmr/data")

# Analyze 2D relaxation data
relaxation2d("expno", "vdlist")

# Analyze heteronuclear NOE data
hetnoe2d(["reference/pdata/1", "saturated/pdata/1"], [false, true])

# Run diffusion analysis
diffusion()

# Run TRACT analysis
tract()
```

## Key Functions

### 2D Analysis with Interactive GUI

- `relaxation2d(experimentfiles, relaxationtimes)`: Analyze T1/T2 relaxation data
- `hetnoe2d(planefilenames, saturation)`: Analyze heteronuclear NOE data

### 1D Analysis with Guided Interface

- `diffusion([coherence])`: Analyze diffusion experiments
- `tract()`: Analyze TRACT experiments for estimating rotational correlation times

### Utility Functions

- `viscosity(solvent, T)`: Calculate solvent viscosity at given temperature
  - Supports `:h2o` and `:d2o`
  - Uses model from Cho et al, J Phys Chem B (1999) 103 1991-1994

## Data Input

The package supports various Bruker data formats:
- Single experiment directories (e.g., ".../11")
- Multiple experiment directories
- Processed data directories (e.g., ".../pdata/1")
- Variable delay lists (for relaxation experiments)

## Usage Example

```julia
# Analyze T2 relaxation data from multiple processed spectra
relaxation2d(
    ["expno/pdata/231",  # First spectrum
     "expno/pdata/232",  # Second spectrum
     "expno/pdata/233"], # Third spectrum
    [0.01, 0.03, 0.05]   # Relaxation delays in seconds
)

# Analyze hetNOE data from reference/saturated pair
hetnoe2d(
    ["reference/pdata/1",
     "saturated/pdata/1"],
    [false, true]
)
```

## Interactive Analysis

The 2D analysis functions provide an interactive GUI for:
- Peak picking
- Fitting
- Visual validation
- Data export

The 1D analysis functions provide guided command-line interfaces for:
- Parameter input
- Region selection
- Data fitting
- Figure generation

## Output

All analysis functions provide:
- Fitted parameters with uncertainties
- Publication-quality figures
- Options to save results

