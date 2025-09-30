# NMRAnalysis.jl

A Julia package for analysis of NMR data, with particular focus on relaxation, diffusion, screening and
protein dynamics experiments.

!!! note
    This package is under active development and its features and syntax may change.

NMRAnalysis.jl provides tools for analysing various types of NMR experiments commonly
used in biomolecular NMR:

- 1D diffusion (`diffusion()`)
- 1D TRACT (`tract()`)
- 1D 19F R1rho relaxation dispersion (`r1rho()`)
- 2D relaxation measurements (T1, T2) (`relaxation2d(experimentfiles, relaxationtimes)`)
- Heteronuclear NOE (`hetnoe2d(planefilenames, saturation)`)

Some utility functions are also provided:

- `viscosity(solvent, T)`: Estimate solution viscosity at given temperature. Solvent can be `:h2o` or `:d2o`.

