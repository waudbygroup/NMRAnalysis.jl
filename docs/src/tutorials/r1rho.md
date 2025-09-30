# Tutorial: R1ρ Analysis

This tutorial walks through the graphical interface for analysing R1ρ relaxation dispersion data. It covers experimental setup, GUI controls, fitting workflow, and interpretation of results.

## Experimental Setup 

### 1. Determine the Ligand Chemical Shift 🔎

Acquire a 1D broadband <sup>19</sup>F NMR experiment to identify the chemical shift of the ligand resonance. This value will serve as the O1 in the 90° pulse calibration experiment. Ensure that the pulse length (P1) is set appropriately in subsequent $R_{1ρ}$ experiments.

### 2. Pulse Program

The pulse program used for on-resonance $R_{1ρ}$ experiments is available on GitHub:

- **[`19F_onresR1p.cw`](https://github.com/chriswaudby/pp/blob/master/19F_onresR1p.cw)**  

This sequence is adapted from Overbeck et al. (J. Magn. Reson. 2020, **74**, 753–766) and includes improved temperature compensation. It supports pseudo-3D acquisition with variable spinlock lengths and powers, read in via `VPLIST` and `VALIST`.

> Make sure to calibrate your hard pulse and spinlock powers before running the sequence. See the calibration section below for details.

### 3. Calibrate Spinlock Powers

NMRAnalysis.jl provides the function `setupR1rhopowers()` for calculating spinlock powers based on a calibrated <sup>19</sup>F 90 degree hard pulse:

- Input the P1 value for the <sup>19</sup>F hard pulse (in μs) and the PLdB1 (<sup>19</sup>F hard pulse power in dB). 
- You may supply a custom list of spinlock strengths (in Hz), or use the default set provided.
- For high spinlock powers, a warning will be issued to verify that the spinlock duration remains within acceptable power limits. 
- The output is a list of calibrated spinlock powers (in dB) that can be copied directly into VALIST in TopSpin. 

![setupR1rhopowers](../assets/setupR1rhopowers.png)

## GUI Overview 

Open a Julia session in the terminal and launch the analysis interface:

```julia
using NMRAnalysis

# prompt to select experiments
r1rho()

# provide path to directory containing experiments
r1rho("example/R1rho")

# provide path to specific experiments
r1rho(["example/R1rho/11", "example/R1rho/12"])
```

 Once loaded, the GUI displays the first spectrum of the dataset.

![R1ρ Analysis Interface](../assets/r1rho-interface.png)

- **Series Toggle**: Switch between measurements at different spin-lock field strengths.
- **Integration Width**: Manually input a value or click **Optimise** to automatically minimize fitting error.
- **Peak Position (ppm)**: Automatically set to the chemical shift of a ligand; manually adjust accordingly if analsying a mixture.
- **Noise Position (ppm)**: Automatically placed away from the peak; adjust if baseline noise is misestimated.
- **Initial Guesses**: Provide starting values for `R2,0`, `Rex`, and `kex` to guide model fitting.
- **Δδ stdev (ppm)**: Accounts for uncertainty in the chemical shift difference between free and bound states. Assumes a normal distribution centered at 0 ppm with a standard deviation of 2 ppm.
- **Output Folder**: Specify a name for your results folder to keep outputs organised.
- **Save Results**: Export fitted parameters and plots to the output folder.

## Analysis Workflow

### 1. Visualise the Spectrum

- The top panel displays the observed spectrum at a given spin-lock field strength ($ν_{SL}$).
- Peak and noise positions are marked and can be adjusted by dragging.
- Set the **Integration Width** to define the region used for peak fitting.

### 2. Fit the Data

- Click **Optimise** to refine the integration width automatically.
- Input initial guesses for model parameters:
  - `R2,0`: Baseline transverse relaxation rate
  - `Rex`: Exchange contribution to relaxation
  - `kex`: Exchange rate constant
- The GUI fits the data and overlays model curves on the plots.

### 3. Interpret the Results

Signal intensities are fit globally as a function of relaxation time and spin-lock field strength:

$$
I(T_{\text{SL}}, \nu_{\text{SL}}) = I_0 \cdot \exp\left(-\left[R_{2,0} + \frac{R_{\text{ex}} \cdot K^2}{K^2 + 4\pi^2 \nu_{\text{SL}}^2}\right] \cdot T_{\text{SL}}\right)
$$

*Adapted from Trott & Palmer (2002), J. Magn. Reson. 154, 157–160.*

Where:

$$
K^2 = k_{\text{ex}}^2 + 4\pi^2 \Delta\nu^2
$$

To assess whether exchange contributes significantly, a null model excluding `Rex` is also fit and compared using an F-test.

#### Dispersion Curve

The GUI plots $R_{1ρ}$ as a function of $ν_{SL}$ using fitted parameters:

$$
R_{1\rho} = R_{2,0} + \frac{R_{\text{ex}} \cdot K^2}{K^2 + 4\pi^2 \nu_{\text{SL}}^2}
$$

This curve is overlaid with $R_{1ρ}$ values obtained from exponential fits at individual spin-lock field strengths, enabling visual comparison of model performance.

#### Chemical Shift Correction

To account for uncertainty in the chemical shift difference ($Δδ$), a particle-based Monte Carlo correction is applied to $K$. For each particle, $k_{\text{off}}$ is calculated as:

$$
k_{\text{off}} \approx k_{\text{ex}} = \sqrt{K^2 - 4\pi^2 \Delta\nu^2}
$$

Samples yielding nonphysical values are excluded, and the final estimate is reported as the mean ± standard deviation of valid particles.

### 4. Output Files

Upon saving results, the GUI generates both raw and fitted data files, along with summary plots. These outputs are organized by analysis type:

#### Dispersion Curve Outputs

These files correspond to the global fit of $R_{1ρ}$ versus spinlock field strength:

- `dispersion-points.csv`: Raw $R_{1ρ}$ values from exponential fits at each spinlock field
- `dispersion-fit.csv`: Fitted $R_{1ρ}$ values based on the global model
- `dispersion.pdf`: Plot of the dispersion curve with overlaid model fit

![dispersion.pdf](../assets/example_results_file.png)

> The CSV files can be loaded into data frames, allowing efficient customization of plots and further exploration of fitted parameters.


#### Peak Integral Outputs

For each spinlock power, the GUI exports:

- `intensities_<spinlock>Hz-points.csv`: Raw peak intensities as a function of spinlock duration
- `intensities_<spinlock>Hz-fit.csv`: Fitted intensities using the relaxation model
- `intensities_<spinlock>Hz.pdf`: Plot of intensity decay curves with fitted overlays

> These files support detailed inspection of signal decay and fitting quality at individual spinlock powers.
 
 #### Summary File: `results.txt`

This file provides a concise summary of the analysis, including:

- Input file paths used in the fitting  
- Peak and noise positions (ppm)  
- Integration width  
- Initial parameter guesses (`I0`, `R2,0`, `Rex`, `kex`)  
- Final fitted values with uncertainties  

**Example:**

![results_file](../assets/example_results_file.png)

> This file is useful for quick reference and record-keeping, especially when comparing fits across multiple datasets or conditions.


