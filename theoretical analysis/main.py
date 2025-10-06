
"""
Brownian Motion Simulations
============================

This project explores Brownian dynamics in one dimension through numerical simulations
of the Langevin equation. Two cases are studied:
1. Free diffusion (no external potential)
2. Harmonic confinement (harmonic oscillator potential)

Author: Lorena Puhl
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from typing import Tuple, List
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
plt.style.use('seaborn-v0_8-poster')
plt.rcParams['figure.figsize'] = (10, 8)


class BrownianSimulation:
    """
    A class to simulate and analyze Brownian motion in 1D.
    
    The Langevin equation is solved numerically:
        dx/dt = -∇U(x) + √(2D) * ξ(t)
    
    where:
        U(x) is the potential energy
        D is the diffusion coefficient
        ξ(t) is white Gaussian noise
    """
    
    def __init__(self, diffusion_coeff: float = 1.0, dt: float = 0.01):
        """
        Initialize the simulation parameters.
        
        Parameters
        ----------
        diffusion_coeff : float
            Diffusion coefficient D (default: 1.0)
        dt : float
            Time step for numerical integration (default: 0.01)
        """
        self.D = diffusion_coeff
        self.dt = dt
        
    def simulate_trajectories(
        self,
        n_simulations: int,
        t_final: float,
        x_initial: float = 0.0,
        force_function=None
    ) -> Tuple[np.ndarray, List[float]]:
        """
        Simulate multiple Brownian trajectories.
        
        Parameters
        ----------
        n_simulations : int
            Number of independent trajectories to simulate
        t_final : float
            Final simulation time
        x_initial : float
            Initial position (default: 0.0)
        force_function : callable or None
            Function computing force F(x) = -dU/dx. If None, free diffusion.
            
        Returns
        -------
        time_grid : np.ndarray
            Array of time points
        final_positions : list
            List of final positions for each trajectory
        """
        # Create time grid
        time_grid = np.arange(0, t_final, self.dt)
        n_steps = len(time_grid)
        
        # Store final positions
        final_positions = []
        
        # Run simulations
        for _ in range(n_simulations):
            # Initialize position array
            x = np.zeros(n_steps)
            x[0] = x_initial
            
            # Evolve according to Langevin equation
            for i in range(n_steps - 1):
                # Generate random noise
                noise = np.random.normal(0, 1)
                
                # Compute force term (if potential present)
                force_term = 0.0 if force_function is None else force_function(x[i])
                
                # Update position: Euler-Maruyama method
                x[i + 1] = x[i] + force_term * self.dt + np.sqrt(2 * self.D * self.dt) * noise
            
            final_positions.append(x[-1])
        
        return time_grid, final_positions
    
    @staticmethod
    def analytical_distribution(x: np.ndarray, t: float, x0: float, D: float) -> np.ndarray:
        """
        Analytical solution for free diffusion probability distribution.
        
        Parameters
        ----------
        x : np.ndarray
            Position values
        t : float
            Time
        x0 : float
            Initial position
        D : float
            Diffusion coefficient
            
        Returns
        -------
        probability : np.ndarray
            Probability density at positions x
        """
        return (1 / np.sqrt(4 * np.pi * D * t)) * np.exp(-((x - x0)**2) / (4 * D * t))
    
    @staticmethod
    def compute_statistics(positions: List[float]) -> Tuple[float, float]:
        """
        Compute mean and variance of position distribution.
        
        Parameters
        ----------
        positions : list
            List of positions
            
        Returns
        -------
        mean : float
            Mean position
        variance : float
            Variance of positions
        """
        positions_array = np.array(positions)
        mean = np.mean(positions_array)
        variance = np.var(positions_array)
        return mean, variance
    
    @staticmethod
    def compute_msd(positions: List[float], reference_position: float = 0.0) -> float:
        """
        Compute Mean Square Displacement.
        
        Parameters
        ----------
        positions : list
            List of final positions
        reference_position : float
            Reference position (default: 0.0)
            
        Returns
        -------
        msd : float
            Mean square displacement
        """
        positions_array = np.array(positions)
        return np.mean((positions_array - reference_position)**2)


    @staticmethod
    def analytical_harmonic_equilibrium(x: np.ndarray, k: float = 1.0, D: float = 1.0) -> np.ndarray:
        """
        Analytical equilibrium distribution for harmonic potential.
        
        For U(x) = 0.5 * k * x², the Boltzmann distribution is:
        P_eq(x) = sqrt(k/(2πD)) * exp(-k*x²/(2D))
        
        Parameters
        ----------
        x : np.ndarray
            Position values
        k : float
            Spring constant (default: 1.0)
        D : float
            Diffusion coefficient (default: 1.0)
            
        Returns
        -------
        probability : np.ndarray
            Equilibrium probability density
        """
        return np.sqrt(k / (2 * np.pi * D)) * np.exp(-k * x**2 / (2 * D))


def plot_distribution(
    positions: List[float],
    time: float,
    analytical_func=None,
    title: str = "Position Distribution",
    x0: float = 0.0,
    D: float = 1.0
):
    """
    Plot histogram of positions and compare with analytical solution.
    
    Parameters
    ----------
    positions : list
        List of final positions
    time : float
        Simulation time
    analytical_func : callable or None
        Function to compute analytical solution
    title : str
        Plot title
    x0 : float
        Initial position
    D : float
        Diffusion coefficient
    """
    positions_array = np.array(positions)
    mean, variance = BrownianSimulation.compute_statistics(positions)
    
    plt.figure(figsize=(10, 8))
    
    # Plot histogram
    n, bins, _ = plt.hist(
        positions_array,
        bins=50,
        density=True,
        alpha=0.7,
        label=f'Simulation\nMean = {mean:.4f}\nVariance = {variance:.4f}'
    )
    
    # Plot analytical solution if provided
    if analytical_func is not None:
        x_sorted = np.sort(positions_array)
        analytical = analytical_func(x_sorted, time, x0, D)
        plt.plot(x_sorted, analytical, 'r.', markersize=3, label='Analytical Solution')
    
    plt.xlabel('Final Position x', fontsize=15)
    plt.ylabel('Probability Density', fontsize=15)
    plt.title(f'{title} (t = {time})', fontsize=16)
    plt.legend(loc='best', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_msd(times: np.ndarray, msd_values: np.ndarray, title: str = "Mean Square Displacement"):
    """
    Plot MSD as a function of time.
    
    Parameters
    ----------
    times : np.ndarray
        Array of time points
    msd_values : np.ndarray
        Array of MSD values
    title : str
        Plot title
    """
    plt.figure(figsize=(10, 8))
    plt.plot(times, msd_values, 'o-', markersize=8, linewidth=2, label='Numerical')
    
    # Add theoretical line for free diffusion (MSD = 2Dt)
    if "Free" in title:
        theoretical = 2 * times  # Assuming D=1
        plt.plot(times, theoretical, '--', linewidth=2, label='Theory: MSD = 2Dt')
    
    plt.xlabel('Time', fontsize=15)
    plt.ylabel('MSD', fontsize=15)
    plt.title(title, fontsize=16)
    plt.legend(loc='best', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def compare_with_equilibrium(
    positions: List[float],
    time: float,
    D: float = 1.0,
    k: float = 1.0
) -> dict:
    """
    Compare simulation results with theoretical equilibrium for harmonic potential.
    
    Parameters
    ----------
    positions : list
        Simulated final positions
    time : float
        Simulation time
    D : float
        Diffusion coefficient
    k : float
        Spring constant
        
    Returns
    -------
    results : dict
        Dictionary containing comparison metrics
    """
    positions_array = np.array(positions)
    
    # Theoretical predictions
    theoretical_variance = D / k
    theoretical_mean = 0.0
    
    # Relaxation time
    tau = 1 / k  # For overdamped dynamics with friction coefficient = 1
    
    # Time-dependent variance (before equilibration)
    variance_t = theoretical_variance * (1 - np.exp(-2 * k * time))
    
    # Empirical statistics
    empirical_mean = np.mean(positions_array)
    empirical_variance = np.var(positions_array)
    
    # Compute relative errors
    mean_error = abs(empirical_mean - theoretical_mean)
    variance_error_eq = abs(empirical_variance - theoretical_variance) / theoretical_variance
    variance_error_t = abs(empirical_variance - variance_t) / variance_t
    
    return {
        'time': time,
        'relaxation_time': tau,
        'empirical_mean': empirical_mean,
        'empirical_variance': empirical_variance,
        'theoretical_mean': theoretical_mean,
        'theoretical_variance_eq': theoretical_variance,
        'theoretical_variance_t': variance_t,
        'mean_error': mean_error,
        'variance_error_equilibrium': variance_error_eq,
        'variance_error_time_dependent': variance_error_t,
        'is_equilibrated': time > 3 * tau
    }


def main():
    """
    Main function to run simulations for both free diffusion and harmonic potential.
    """
    # Simulation parameters
    D = 1.0
    dt = 0.01
    n_simulations = 10000
    time_points = np.array([0.5, 1.0, 3.0, 5.0, 10.0])
    x0 = 0.0
    
    # Initialize simulator
    simulator = BrownianSimulation(diffusion_coeff=D, dt=dt)
    
    print("=" * 60)
    print("BROWNIAN MOTION SIMULATIONS")
    print("=" * 60)
    
    # ========== FREE DIFFUSION ==========
    print("\n1. FREE DIFFUSION (No External Potential)")
    print("-" * 60)
    
    msd_free = []
    
    for t_final in time_points:
        print(f"\nSimulating {n_simulations} trajectories until t = {t_final}...")
        
        _, final_positions = simulator.simulate_trajectories(
            n_simulations=n_simulations,
            t_final=t_final,
            x_initial=x0,
            force_function=None  # Free diffusion
        )
        
        # Compute MSD
        msd = simulator.compute_msd(final_positions, reference_position=x0)
        msd_free.append(msd)
        
        mean, variance = simulator.compute_statistics(final_positions)
        print(f"  Mean = {mean:.4f}, Variance = {variance:.4f}, MSD = {msd:.4f}")
        
        # Plot distribution
        plot_distribution(
            positions=final_positions,
            time=t_final,
            analytical_func=simulator.analytical_distribution,
            title="Free Diffusion - Position Distribution",
            x0=x0,
            D=D
        )
    
    # Plot MSD vs time for free diffusion
    plot_msd(time_points, np.array(msd_free), title="Free Diffusion - MSD vs Time")
    
    # ========== HARMONIC POTENTIAL ==========
    print("\n2. HARMONIC POTENTIAL (U = 0.5 * k * x²)")
    print("-" * 60)
    
    # Define force for harmonic potential: F = -kx (with k=1)
    def harmonic_force(x):
        return -x  # k = 1, F = -dU/dx = -kx
    
    msd_harmonic = []
    
    for t_final in time_points:
        print(f"\nSimulating {n_simulations} trajectories until t = {t_final}...")
        
        _, final_positions = simulator.simulate_trajectories(
            n_simulations=n_simulations,
            t_final=t_final,
            x_initial=x0,
            force_function=harmonic_force
        )
        
        # Compute MSD
        msd = simulator.compute_msd(final_positions, reference_position=x0)
        msd_harmonic.append(msd)
        
        mean, variance = simulator.compute_statistics(final_positions)
        print(f"  Mean = {mean:.4f}, Variance = {variance:.4f}, MSD = {msd:.4f}")
        
        # Plot distribution
        plot_distribution(
            positions=final_positions,
            time=t_final,
            analytical_func=lambda x, t, x0, D: simulator.analytical_harmonic_equilibrium(x, k=1.0, D=D),
            title="Harmonic Potential - Position Distribution",
            x0=x0,
            D=D
        )
    
    # Plot MSD vs time for harmonic potential
    plot_msd(time_points, np.array(msd_harmonic), title="Harmonic Potential - MSD vs Time")

    # Compare simulation results with theoretical equilibrium for harmonic potential
    comparison = compare_with_equilibrium(final_positions, t_final, D=D, k=1.0)

    print(f"  Mean = {comparison['empirical_mean']:.4f} (theory: {comparison['theoretical_mean']:.4f})")
    print(f"  Variance = {comparison['empirical_variance']:.4f}")
    print(f"  Theory (equilibrium): {comparison['theoretical_variance_eq']:.4f}")
    print(f"  Theory (time-dependent): {comparison['theoretical_variance_t']:.4f}")
    print(f"  Equilibrated: {comparison['is_equilibrated']}")
    print(f"  Relative error: {comparison['variance_error_time_dependent']*100:.2f}%")
    
    print("\n" + "=" * 60)
    print("SIMULATIONS COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Run simulations
    main()




