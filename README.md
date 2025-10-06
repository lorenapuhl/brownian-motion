# Brownian Motion Analysis Repository

A collection of computational physics projects exploring Brownian dynamics through both theoretical simulations and experimental data analysis.

## Introduction

This repository contains two complementary approaches to studying Brownian motion:

1. **Theoretical Analysis**: Monte Carlo simulations of the Langevin equation for free diffusion and harmonic confinement
2. **Experimental Analysis**: Statistical analysis of real experimental trajectories from particle tracking experiments

Both projects investigate fundamental questions in statistical mechanics, stochastic processes, and non-equilibrium physics.

## Repository Structure

```
brownian-motion-analysis/
│
├── README.md                    # This file - general project overview
│
├── theoretical_analysis/
│   ├── README.md                # Detailed documentation for simulations
│   ├── main.py                  # Main simulation script
│   └── results/
│       ├── free_diffusion.png           # Position distributions at various times
│       ├── free_diffusion_msd.png       # MSD vs time for free diffusion
│       ├── harmonic_potential.png       # Position distributions in harmonic trap
│       └── harmonic_potential_msd.png   # MSD vs time for confined motion
│
└── experimental_analysis/
    ├── README.md                # Documentation for experimental data analysis
    ├── main.py                  # Data processing and analysis script
    └── results/
        ├── free_diffusion.png           # Experimental trajectory analysis
        ├── free_diffusion_msd.png       # Measured MSD from experiments
        ├── harmonic_potential.png       # Confined particle trajectories
        └── harmonic_potential_msd.png   # MSD in experimental harmonic trap
```

## Projects Overview

### Theoretical Analysis
Numerical simulations solving the Langevin equation:
- Free Brownian diffusion (no external forces)
- Harmonic confinement (quadratic potential)
- Validation against analytical solutions
- Mean Square Displacement analysis

**Key Results**: Linear MSD growth for free diffusion, MSD saturation for confinement

### Experimental Analysis
Statistical analysis of particle tracking data:
- Real experimental trajectories
- Comparison with theoretical predictions
- Parameter estimation (diffusion coefficient, trap stiffness)
- Verification of stochastic theory

**Key Results**: Quantitative agreement between experiment and theory

## Quick Start

Each subdirectory contains its own detailed README with:
- Theoretical background
- Usage instructions
- Results and plots
- Detailed analysis

Navigate to `theoretical_analysis/` or `experimental_analysis/` to get started.

## Requirements

```bash
pip install numpy matplotlib scipy
```

Python 3.7+

## Author

Lorena Puhl

lorena.puhl@protonmail.com

## License

MIT License
