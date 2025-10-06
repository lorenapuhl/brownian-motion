# Brownian Motion Analysis

A Python script for analyzing experimental Brownian motion data to calculate the Boltzmann constant and diffusion coefficient.

## Theory

Brownian motion is the random movement of particles in a fluid caused by collisions with molecules. The Boltzmann constant (k<sub>B</sub>) links particle motion to temperature, while the diffusion coefficient (D) describes how quickly particles spread over time. Together, these concepts help connect microscopic particle dynamics to macroscopic behavior.

<p align="center">
  <img src="https://github.com/user-attachments/assets/f222eccb-5129-413c-b1a3-036e05672c83" alt="Displacement Histogram" width="800"/>
     <br>
  <em>Figure 1: Particle trajectories determined by brownian motion. (Figure from https://towardsdatascience.com/stochastic-processes-simulation-brownian-motion-the-basics-c1d71585d9f9/) </em>
</p>
 
According to the **Einstein–Smoluchowski relation**, the diffusion coefficient $\( D \)$ is given by:

<p align="center">
$$
D = \frac{k_B T}{6\pi \eta r}
$$


where  
- $k_B$ — Boltzmann constant  
- $T$ — absolute temperature  
- $\eta$ — dynamic viscosity of the fluid  
- $r$ — radius of the particle  

From this, the **Boltzmann constant** can be determined as:

$$
k_B = \frac{6 \pi \eta r D}{T}
$$
</p>

---

## Overview

This script processes particle trajectory data from Brownian motion experiments and performs statistical analysis to determine fundamental physical constants. The analysis is based on the Einstein relation for diffusion.

## Repository Structure

```
├── README.md                    # This file - project overview
├── main.py                      # python script for analysis
├── measurements.dat             # experimental data
├── results/                     # analysis results
│   ├── Graph_Brown.png          # Visualising particle displacements
│   ├── Brown_histogram.png      # Histogram of particle displacements
│   ├── Brown_cumulative.png     # Plot cumulative displacement over time

```


## Requirements

```bash
numpy
matplotlib
scipy
```

Install dependencies:
```bash
pip install numpy matplotlib scipy
```

## Input Data Format

The script expects a tab/space-separated data file with the following columns:
- Column 1: Time (seconds)
- Column 2: X position (μm)
- Column 3: Y position (μm)

Note: Decimal values should use commas (e.g., "1,234") - the script automatically converts these to floats.

## Configuration

Before running, adjust the following parameters in the script:

```python
DATA_FILE = 'measurements.dat'  # Path to your data file
TEMPERATURE = 23.4 + 273.15     # Temperature in Kelvin
VISCOSITY = 9.2e-4              # Viscosity in Pa·s
PARTICLE_RADIUS = 0.5 * 755e-9  # Particle radius in meters
```

## Usage

```bash
python main.py
```

## Output

The script generates:

1. **Console output:**
   - Data loading summary
   - Calculated Boltzmann constant with uncertainty
   - Calculated diffusion coefficient with uncertainty
   - Comparison with literature values

2. **Three PDF figures:**
   - `Graph_Brown.pdf`: Particle trajectory plot (x vs y position)
   - `Brown_histogram.pdf`: Histogram of displacements with Gaussian fit
   - `Brown_cumulative.pdf`: Cumulative squared displacement vs time

## Analysis Methods

The script performs the following analysis steps:
1. **Trajectory Visualization**: Plots the 2D path of the particle.  
2. **Displacement Calculation**: Computes step-by-step position changes.  
3. **Mean Squared Displacement**: Calculates $\langle r^2 \rangle$ and time intervals.  
4. **Boltzmann Constant**: Determines $k_B$ using:  
   $k_B = \frac{6 \pi \eta a \langle r^2 \rangle}{4 T \langle \Delta t \rangle}$
5. **Diffusion Coefficient**: Calculates  
   $D = \frac{k_B T}{6 \pi \eta a}$
6. **Statistical Validation**: Creates histogram and cumulative displacement plots.

---

## Results

### Calculated Physical Constants

The analysis of the experimental data yielded:

- **Boltzmann constant**: (1.7680e-23 ± 4.3588e-24) J/K
- **Literature value**: 1.381 × 10⁻²³ J/K
- **Diffusion coefficient**: (8.0091e-13 ± 2.0019e-13) m²/s

### Generated Figures

#### 1. Particle Trajectory

<p align="center">
  <img src="https://github.com/user-attachments/assets/1a5873fa-1694-4d3c-a365-09066599bf7f" alt="Brownian Motion Trajectory" width="800"/>
</p>

The trajectory shows the random walk of the particle due to collisions with surrounding molecules.

#### 2. Displacement Distribution
<p align="center">
  <img src="https://github.com/user-attachments/assets/b38b9106-8ec0-406a-b3d7-1045b68ec79a" alt="Displacement Histogram" width="800"/>
</p>


The histogram demonstrates that displacements follow a Gaussian distribution, as predicted by the theory of Brownian motion. The fitted Gaussian parameters are:
- μ = -0.06 μm
- σ = 1.27 μm

#### 3. Cumulative Squared Displacement
<p align="center">
  <img src="https://github.com/user-attachments/assets/e0881dc2-794f-49f6-9a74-4b620c1d7428" width="800"/>
</p>
The linear relationship between cumulative squared displacement and time confirms Einstein's diffusion relation. The slope of the fitted line provides an independent measure of the diffusion coefficient.

___
## Author

Lorena Puhl

lorena.puhl@protonmail.com
