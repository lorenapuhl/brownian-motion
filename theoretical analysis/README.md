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

<p align="center">

  <img src="https://github.com/user-attachments/assets/e07d783d-972c-4b24-8c80-cabed14d2889" alt="Image 1" width="700"/>
  <br>
  <em>Figure 1: Brownian motion of a particle under free diffusion. (Image from https://www.quantstart.com/articles/brownian-motion-simulation-with-python/). </em>
</p>

2. **Harmonic confinement** - particle trapped in a quadratic potential well (i.e. Optical tweezers experiments, ion traps, polymer dynamics near equilibrium)

The simulations verify key theoretical predictions from statistical mechanics and stochastic processes.

## Theoretical Background

### The Langevin Equation

Brownian motion is described by the overdamped Langevin equation:

$$
\frac{dx}{dt} = F(x) + \sqrt{2D} \xi(t)
$$

where:
- $x(t)$ is the particle position at time $t$
- $F(x) = -\frac{dU}{dx}$ is the deterministic force derived from the potential $U(x)$
- $D$ is the diffusion coefficient, related to temperature and friction. It quantifies how rapidly particles spread out over time due to random thermal motion in a medium.
- $\xi(t)$ is white Gaussian noise with $\langle \xi(t) \rangle = 0$ and $\langle \xi(t)\xi(t') \rangle = \delta(t - t')$

### Case 1: Free Diffusion (U = 0)

For a free particle with no external forces (F = 0):

### Case 1: Free Diffusion

**Analytical solution:**

$$
P(x,t) = \frac{1}{\sqrt{4 \pi D t}} \exp\left[-\frac{(x - x_0)^2}{4 D t}\right]
$$

**Key predictions:**
- Position distribution is Gaussian with variance $\( \sigma^2 = 2Dt \)$
- Mean Square Displacement: **$MSD = 2Dt$** (linear growth)
- Standard deviation grows as $\sqrt{t}$

---

### Case 2: Harmonic Potential $U(x) = \tfrac{1}{2} k x^2$

For a particle in a harmonic trap with force $F = -kx$ :

**Equilibrium distribution (Boltzmann):**

$$
P_{\text{eq}}(x) = \sqrt{\frac{k}{2\pi D}} \exp\left(-\frac{kx^2}{2D}\right)
$$

**Key predictions:**
- Equilibrium variance: $\sigma^2_{\text{eq}} = \frac{D}{k}$
- Relaxation time: $\tau = \frac{1}{k}$
- Time-dependent variance: $\sigma^2(t) = \frac{D}{k}\left[1 - e^{-2kt}\right]$
- **MSD saturates** to constant value $\sim \frac{2D}{k}$ at long times  
- At short times ($t \ll \tau$): diffusive behavior  
- At long times ($t \gg \tau$): confined behavior  


---

### Numerical Method

The Langevin equation is integrated using the **Euler–Maruyama scheme**:

$$
x_{i+1} = x_i + F(x_i)\,\Delta t + \sqrt{2D\Delta t}\,W_i
$$


where $W_i \sim N(0,1)$ are independent Gaussian random variables.

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
python main.py
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
from main import BrownianSimulation

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

<p align="center">
  <img src="https://github.com/user-attachments/assets/4fefcb28-b949-467e-8733-dd5cde625aa0"  alt="Image 2" width="45%"/>
  <img src="https://github.com/user-attachments/assets/8fd68781-6676-4012-9f43-a4992600b83d"alt="Image 1" width="45%"/>
</p>

<p align="center">
  <em>Figure 3: Position distributions at t = 0.5 (left) and t= 10.0. (right). Histograms show simulation results, red dots show analytical Gaussian solution.</em>
</p>



<p align="center">
  <img src="https://github.com/user-attachments/assets/8470cc97-f195-4ee5-8fea-6db4c280fe87" width="700"/>
  </p>

<p align="center">
  <em>Figure 4: Mean Square Displacement grows linearly with time (MSD = 2Dt), confirming diffusive behavior.</em>
</p>

**Key observations:**
- Distributions remain Gaussian at all times
- Width increases as $\sqrt{t}$
- Excellent agreement with analytical predictions
- MSD slope ≈ 2D (for $D=1$, slope ≈ 2)

### Harmonic Potential

<p align="center">
  <img src="https://github.com/user-attachments/assets/39845f42-b219-47b6-8638-b46cc2eabb59" width="45%"/>
  <img src="https://github.com/user-attachments/assets/f1533d3b-43d6-42ca-9176-6599ea0d4375" width="45%"/>
    </p>

<p align="center">
  <em>Figure 5: Position distributions at t = 0.5 (left) and t= 10.0. (right).</em>
</p>


<p align="center">
  <img src="https://github.com/user-attachments/assets/46beabf4-1d00-4439-a0ab-edac66d53808" width="700"/>
  </p>

<p align="center">
  <em>Figure 6: Position distributions approach equilibrium Boltzmann distribution (Gaussian with  σ^2 = D / k). MSD saturates to ~1 at long times, showing confinement effect.</em>
</p>

**Key observations:**
- Initial diffusive regime (t < τ)
- Convergence to equilibrium distribution at t ~ 3τ
- MSD saturation confirming harmonic confinement
- Equilibrium variance ≈ 1 (for D=1, k=1)

#### Statistical Comparison


Simulating 10000 trajectories until t = 10.0...

| Time | Mean   | Variance (Sim) | Variance (Theory) | Relative Error |
|------|--------|----------------|-------------------|----------------|
| 10.0 | 0.0076 | 0.9828         | 1.0000            | 1.72%          |

---

## Theoretical Insights

### 1. Universality of Diffusion
The simulations demonstrate that free Brownian motion leads to a universal diffusive spreading with MSD ∝ t, regardless of microscopic details. This behavior is fundamental to many physical, chemical, and biological processes.

### 2. Role of Confining Potentials
The harmonic potential case illustrates how external forces fundamentally alter transport properties:
- **Transition from transport to confinement**: At short times, particles diffuse freely. At long times, the restoring force dominates, leading to localization.
- **Equilibrium fluctuations**: Even at equilibrium, thermal noise causes persistent fluctuations with amplitude set by the balance between D (driving) and k (restoring).


### 3. Fluctuation-Dissipation Connection
The equilibrium variance $\sigma^2_{eq} = D/k$ embodies the *fluctuation-dissipation theorem*: stronger dissipation (larger k) leads to smaller equilibrium fluctuations, with thermal energy (D) driving the fluctuations.

### 4. Numerical Method Validation
The agreement between simulations and analytical solutions validates:
- The Euler-Maruyama numerical scheme for stochastic differential equations
- Convergence properties with time step Δt = 0.01
- Statistical accuracy with 10,000 trajectories

---
## Author
Lorena Puhl

lorena.puhl@protonmail.com
