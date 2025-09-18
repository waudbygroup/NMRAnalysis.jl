# NMRAnalysis.jl

A Julia package for analysis of NMR data, with particular focus on relaxation, diffusion, screening and
protein dynamics experiments.

!!! note
    This package is under active development and its features and syntax may change.

## Overview

NMRAnalysis.jl provides tools for analysing various types of NMR experiments commonly
used in biomolecular NMR:

- Relaxation measurements (T1, T2)
- Heteronuclear NOE
- Diffusion
- TRACT (TROSY/anti-TROSY Cross-Correlated Relaxation)
- 19F R1rho relaxation dispersion


## Key Functions

### 1D Analysis with Guided Interface

- `diffusion()`: Analyze diffusion experiments
- `tract()`: Analyze TRACT experiments for estimating rotational correlation times

### 1D R1rho relaxation dispersion

- `r1rho()`: Analyse 1D (on-resonance) R1rho relaxation dispersion experiments with interactive interface.

### 2D Analysis with Interactive GUI

- `relaxation2d(experimentfiles, relaxationtimes)`: Analyze T1/T2 relaxation data
- `hetnoe2d(planefilenames, saturation)`: Analyze heteronuclear NOE data

### Utility Functions

- `viscosity(solvent, T)`: Estimate solvent viscosity at given temperature
  - Supports `:h2o` and `:d2o`
  - Uses model from Cho et al, J Phys Chem B (1999) 103 1991-1994

