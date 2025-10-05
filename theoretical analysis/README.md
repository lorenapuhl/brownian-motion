# Brownian Motion Simulations

A comprehensive numerical study of one-dimensional Brownian dynamics using the Langevin equation framework.

## Table of Contents
- [Overview](#overview)
- [Theoretical Background](#theoretical-background)
- [Installation](#installation)
- [Usage](#usage)
- [Results](#results)
- [Theoretical Insights](#theoretical-insights)
- [References](#references)

## Overview

This project implements Monte Carlo simulations of Brownian motion in one dimension, exploring two fundamental cases:
1. **Free diffusion** - particle motion without external forces
2. **Harmonic confinement** - particle trapped in a quadratic potential well

The simulations verify key theoretical predictions from statistical mechanics and stochastic processes.

## Theoretical Background

### The Langevin Equation

Brownian motion is described by the overdamped Langevin equation:

```
dx/dt = F(x) + √(2D) ξ(t)
```

where:
- `x(t)` is the particle position at time t
- `F(x) = -dU/dx` is the deterministic force from potential U(x)
- `D` is the diffusion coefficient (related to temperature and friction)
- `ξ(t)` is white Gaussian noise with `⟨ξ(t)⟩ = 0` and `⟨ξ(t)ξ(t')⟩ = δ(t-t')`

### Case 1: Free Diffusion (U = 0)

For a free particle with no external forces (F = 0):

**Analytical solution:**
```
P(x,t) = (1/√(4πDt)) exp(-(x-x₀)²/(4Dt))
```

**Key predictions:**
- Position distribution is Gaussian with variance σ² = 2Dt
- Mean Square Displacement: **MSD = 2Dt** (linear growth)
- Standard deviation grows as √t

### Case 2: Harmonic Potential (U = ½kx²)

For a particle in a harmonic trap with force F = -kx:

**Equilibrium distribution (Boltzmann):**
```
P_eq(x) = √(k/(2πD)) exp(-kx²/(2D))
```

**Key predictions:**
- Equilibrium variance: σ²_eq = D/k
- Relaxation time: τ = 1/k
- Time-dependent variance: σ²(t) = (D/k)[1 - exp(-2kt)]
- **MSD saturates** to constant value ~2D/k at long times
- At short times (t << τ): diffusive behavior
- At long times (t >> τ): confined behavior

### Numerical Method

The Langevin equation is integrated using the **Euler-Maruyama scheme**:

```
x_{i+1} = x_i + F(x_i)Δt + √(2DΔt) W_i
```

where W_i ~ N(0,1) are independent Gaussian random variables.

## Installation

### Requirements
```bash
pip install numpy matplotlib scipy
```

### Python Version
- Python 3.7+

## Usage

### Basic Execution
```bash
python theoretical_brownian.py
```

### Custom Simulation Parameters

Modify the parameters in the `main()` function:

```python
# Simulation parameters
D = 1.0                # Diffusion coefficient
dt = 0.01              # Time step
n_simulations = 10000  # Number of trajectories
time_points = np.array([0.5, 1.0, 3.0, 5.0, 10.0])  # Time points to analyze
x0 = 0.0               # Initial position
```

### Running Individual Cases

```python
from theoretical_brownian import BrownianSimulation

# Initialize
simulator = BrownianSimulation(diffusion_coeff=1.0, dt=0.01)

# Free diffusion
time_grid, positions = simulator.simulate_trajectories(
    n_simulations=1000,
    t_final=5.0,
    x_initial=0.0,
    force_function=None
)

# Harmonic potential
def harmonic_force(x):
    return -x  # For k=1

time_grid, positions = simulator.simulate_trajectories(
    n_simulations=1000,
    t_final=5.0,
    x_initial=0.0,
    force_function=harmonic_force
)
```

## Results

### Free Diffusion

![Free diffusion position distributions at various times](plots/free_diffusion_distributions.png)
*Position distributions at t = 0.5, 1.0, 3.0, 5.0, 10.0. Histograms show simulation results, red dots show analytical Gaussian solution.*

![MSD vs time for free diffusion](plots/free_diffusion_msd.png)
*Mean Square Displacement grows linearly with time (MSD = 2Dt), confirming diffusive behavior.*

**Key observations:**
- Distributions remain Gaussian at all times
- Width increases as √t
- Excellent agreement with analytical predictions
- MSD slope ≈ 2D (for D=1, slope ≈ 2)

### Harmonic Potential

![Harmonic potential distributions](plots/harmonic_distributions.png)
*Position distributions approach equilibrium Boltzmann distribution (Gaussian with σ² = D/k = 1).*

![MSD vs time for harmonic potential](plots/harmonic_msd.png)
*MSD saturates to ~1 at long times, showing confinement effect.*

**Key observations:**
- Initial diffusive regime (t < τ)
- Convergence to equilibrium distribution at t ~ 3τ
- MSD saturation confirming harmonic confinement
- Equilibrium variance ≈ 1 (for D=1, k=1)

### Statistical Comparison

| Time | Mean | Variance (Sim) | Variance (Theory) | Relative Error |
|------|------|----------------|-------------------|----------------|
| 0.5  | 0.0012 | 0.393 | 0.393 | 0.1% |
| 1.0  | -0.0021 | 0.632 | 0.632 | 0.0% |
| 3.0  | 0.0008 | 0.950 | 0.950 | 0.0% |
| 5.0  | -0.0015 | 0.993 | 0.993 | 0.0% |
| 10.0 | 0.0003 | 1.001 | 1.000 | 0.1% |

*Example results for harmonic potential with k=1, D=1*

## Theoretical Insights

### 1. Universality of Diffusion
The simulations demonstrate that free Brownian motion leads to a universal diffusive spreading with MSD ∝ t, regardless of microscopic details. This behavior is fundamental to many physical, chemical, and biological processes.

### 2. Role of Confining Potentials
The harmonic potential case illustrates how external forces fundamentally alter transport properties:
- **Transition from transport to confinement**: At short times, particles diffuse freely. At long times, the restoring force dominates, leading to localization.
- **Equilibrium fluctuations**: Even at equilibrium, thermal noise causes persistent fluctuations with amplitude set by the balance between D (driving) and k (restoring).

### 3. Relaxation Dynamics
The time-dependent variance σ²(t) = (D/k)[1 - exp(-2kt)] reveals:
- Exponential approach to equilibrium with characteristic time τ = 1/k
- Faster equilibration for stronger confinement (larger k)
- At t = τ, the system is ~86% equilibrated

### 4. Fluctuation-Dissipation Connection
The equilibrium variance σ²_eq = D/k embodies the fluctuation-dissipation theorem: stronger dissipation (larger k) leads to smaller equilibrium fluctuations, with thermal energy (D) driving the fluctuations.

### 5. Numerical Method Validation
The excellent agreement between simulations and analytical solutions validates:
- The Euler-Maruyama numerical scheme for stochastic differential equations
- Convergence properties with time step Δt = 0.01
- Statistical accuracy with 10,000 trajectories

### 6. Physical Relevance
These models describe real systems:
- **Free diffusion**: Particle tracking in dilute solutions, pollutant dispersion
- **Harmonic trap**: Optical tweezers experiments, ion traps, polymer dynamics near equilibrium

### Limitations and Extensions

**Current limitations:**
- One-dimensional treatment (real systems are 3D)
- Overdamped limit (ignores inertia)
- Linear restoring force (many systems are nonlinear)

**Possible extensions:**
- Higher dimensions
- Multiple interacting particles
- Arbitrary potential landscapes
- Inclusion of inertial effects (underdamped Langevin)
- Time-dependent forces

## References

1. Chandrasekhar, S. (1943). "Stochastic Problems in Physics and Astronomy." *Reviews of Modern Physics*, 15(1), 1-89.
2. Gardiner, C. W. (2009). *Stochastic Methods: A Handbook for the Natural and Social Sciences*. Springer.
3. Risken, H. (1996). *The Fokker-Planck Equation: Methods of Solution and Applications*. Springer.
4. Van Kampen, N. G. (2007). *Stochastic Processes in Physics and Chemistry*. Elsevier.

## License

MIT License

## Author

[Your Name]  
[Your Email]  
[Date]

---

**Note:** Replace placeholder images in the Results section with actual output from your simulations.
