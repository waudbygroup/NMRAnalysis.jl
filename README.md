# NMRAnalysis

[![Build Status](https://github.com/waudbygroup/NMRAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/waudbygroup/NMRAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://waudbygroup.github.io/NMRTools.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://waudbygroup.github.io/NMRTools.jl/dev)
[![CI](https://github.com/waudbygroup/NMRTools.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/waudbygroup/NMRTools.jl/actions/workflows/Runtests.yml)
[![Codecov](https://codecov.io/gh/waudbygroup/NMRTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/waudbygroup/NMRTools.jl)
[![DOI](https://zenodo.org/badge/251587402.svg)](https://zenodo.org/badge/latestdoi/251587402) -->

## Outline

* pseudo-2D analysis
  - tasks
    - diffusion (difframp)
    - relaxation (T1, T2) (vdlist, vclist)
    - kinetics
  - break it down:
    - load data
    - assigning values to the indirect dimension
    - integrating slices
    - fitting to a model (not kinetics)
    - plotting slice(s) and fit(s)
    - reporting