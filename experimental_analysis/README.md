# Brownian Motion Analysis

A Python script for analyzing experimental Brownian motion data to calculate the Boltzmann constant and diffusion coefficient.

Brownian motion is the random movement of particles in a fluid caused by collisions with molecules. The Boltzmann constant (k<sub>B</sub>) links particle motion to temperature, while the diffusion coefficient (D) describes how quickly particles spread over time. Together, these concepts help connect microscopic particle dynamics to macroscopic behavior.

---

## Overview

This script processes particle trajectory data from Brownian motion experiments and performs statistical analysis to determine fundamental physical constants. The analysis is based on the Einstein relation for diffusion.

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

1. **Trajectory Visualization**: Plots the 2D path of the particle
2. **Displacement Calculation**: Computes step-by-step position changes
3. **Mean Squared Displacement**: Calculates ⟨r²⟩ and time intervals
4. **Boltzmann Constant**: Determines k_B using: k_B = (6πηa⟨r²⟩) / (4T⟨Δt⟩)
5. **Diffusion Coefficient**: Calculates D = k_B·T / (6πηa)
6. **Statistical Validation**: Creates histogram and cumulative displacement plots

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
  <img src="https://github.com/lorenapuhl/brownian-motion/blob/main/experimental_analysis/Graph_Brown.pdf" alt="Brownian Motion Trajectory" width="800"/>
</p>

The trajectory shows the random walk of the particle due to collisions with surrounding molecules.

#### 2. Displacement Distribution
<p align="center">
  <img src="https://github.com/lorenapuhl/brownian-motion/blob/main/experimental_analysis/Brown_histogram.pdf" alt="Displacement Histogram" width="800"/>
</p>

The histogram demonstrates that displacements follow a Gaussian distribution, as predicted by the theory of Brownian motion. The fitted Gaussian parameters are:
- μ = -0.06 μm
- σ = 1.27 μm

#### 3. Cumulative Squared Displacement
<p align="center">
  <img src="https://github.com/lorenapuhl/brownian-motion/blob/main/experimental_analysis/Brown_cumulative.pdf" alt="Cumulative Displacement" width="800"/>
</p>
The linear relationship between cumulative squared displacement and time confirms Einstein's diffusion relation. The slope of the fitted line provides an independent measure of the diffusion coefficient.

## Author

Lorena Puhl
