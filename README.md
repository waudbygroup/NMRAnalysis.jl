# NMRAnalysis

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://waudbygroup.github.io/NMRAnalysis.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://waudbygroup.github.io/NMRAnalysis.jl/dev)
[![CI](https://github.com/waudbygroup/NMRAnalysis.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/waudbygroup/NMRAnalysis.jl/actions/workflows/Runtests.yml)
[![Codecov](https://codecov.io/gh/waudbygroup/NMRAnalysis.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/waudbygroup/NMRAnalysis.jl)

NMRAnalysis.jl is a library for analysis of NMR experiments such as diffusion and relaxation.

> **NOTE**: This package is under active development and it may be a while before its features and syntax stabilises.

---
## Getting Started with the R1rho GUI

To use the R1rho GUI, follow these steps:

1. **Load the NMRAnalysis package**:
    ```julia
    using NMRAnalysis
    ```

2. **Launch the R1rho GUI**:
    ```julia
    R1rho("path/to/your/experiment/folder")
    ```

3. **Select your experiment folder**:
    - Choose a numbered experiment folder or parent directory.

4. **Navigate and select experiments**:
    - Use `Enter` to select experiments.
    - Press `d` when done.

---

### Working with the GUI

Once the GUI opens, you can set the peak position and adjust fitting parameters to extract kinetic values such as `kex`.

#### GUI Features

- **Series Toggle**: Navigate between measurements at different spin-lock fields.
- **Optimise**: Automatically adjusts the integration width to minimize error.
- **Initial Guesses**: Adjust values for `R2,0`, `Rex`, and `kex` as needed.
- **Δδ stdev (ppm)**: Corrects for the unknown contribution of Δδ.

---

### Saving Your Work

- Specify the **Output folder**.
- Click **Save results** to export your analysis.
