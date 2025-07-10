# R1ρ Analysis

The `R1rho` module in NMRAnalysis.jl provides a graphical interface for the analysis of one-dimensional R1ρ relaxation dispersion experiments. This GUI allows you to load, visualize, and fit R1ρ data interactively.

## Launching the R1ρ GUI

You can launch the R1ρ analysis GUI in several ways, depending on your workflow and data organization:

### 1. Launch with a Selection Dialog

If you call `r1rho()` with no arguments, a dialog will appear allowing you to select a directory containing your NMR experiments.

```julia
using NMRAnalysis
r1rho()
```

### 2. Launch with a Starting Folder

You can provide a starting folder as an argument. The program will display a list of available NMR experiments in the terminal for you to select.

```julia
r1rho("example/R1rho")
```

![R1ρ Experiment Selection](assets/r1rho-selection.png)

### 3. Launch with a List of Specific Input Spectra

You can also provide a list of specific experiment folders or files for direct analysis.

```julia
r1rho(["example/R1rho/11", "example/R1rho/12"])
```

![R1ρ Analysis Interface](assets/r1rho-interface.png)

## Adjusting the Display Size: `scalefactor`

The GUI display size can be adjusted using the optional `scalefactor` keyword argument. By default, this is set to `:automatic`, which uses a value of `2` for high-resolution (HiDPI/Retina) displays and `1` for standard displays. You can override this by specifying a numeric value:

```julia
r1rho(scalefactor=1.5)  # Scale the display by 1.5x
```

!!! tip
    If the interface is too crowded for your screen size, try running with a lower scale factor.

## Example Usage

```julia
# Launch with dialog
r1rho()

# Launch with a specific folder
r1rho("data/R1rho_experiments")

# Launch with a list of experiments and custom scale
r1rho(["data/R1rho/11", "data/R1rho/12"], scalefactor=1.5)
```

## Citation

*Placeholder for the associated publication. Please insert citation details here.*
