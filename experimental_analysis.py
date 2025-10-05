"""
Brownian Motion Experimental Analysis
======================================
This script analyzes experimental data of Brownian motion to calculate
the Boltzmann constant and diffusion coefficient.

Author: [Your name]
Date: [Date]
"""

# %% Import required libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.optimize import curve_fit
from scipy.stats import norm


# %% Configuration
# File path to the measurement data
DATA_FILE = 'measurements.dat'

# Physical constants and experimental parameters
TEMPERATURE = 23.4 + 273.15  # Temperature in Kelvin (converted from Celsius)
TEMP_UNCERTAINTY = 0.1  # Temperature uncertainty in K
VISCOSITY = 9.2e-4  # Viscosity in Pa·s
VISCOSITY_UNCERTAINTY = 0.1e-4  # Viscosity uncertainty
PARTICLE_RADIUS = 0.5 * 755e-9  # Particle radius in meters (755 nm diameter)
RADIUS_UNCERTAINTY = 0.5 * 30e-9  # Radius uncertainty

# Output directory for saving figures
OUTPUT_DIR = "C:/Users/lorena/Desktop/GitHub"


# %% Helper function to convert comma decimals to float
def comma_to_float(value_str):
    """
    Convert comma-separated decimal strings to float.
    
    Parameters:
    -----------
    value_str : bytes
        String with comma as decimal separator
        
    Returns:
    --------
    float
        Converted numerical value
    """
    return float(value_str.decode("utf-8").replace(',', '.'))


# %% Section 1: Load and visualize trajectory data
print("="*50)
print("Section 1: Loading and visualizing trajectory data")
print("="*50)

# Load data: time (column 1), x-position (column 2), y-position (column 3)
time, x_pos, y_pos = np.loadtxt(
    DATA_FILE, 
    skiprows=1,  # Skip header row
    usecols=(1, 2, 3),  # Use columns 1, 2, 3
    converters={1: comma_to_float, 2: comma_to_float, 3: comma_to_float},
    unpack=True
)

print(f"Loaded {len(time)} data points")
print(f"Time range: {time[0]:.2f} to {time[-1]:.2f} seconds")

# Plot particle trajectory
plt.figure(figsize=(11.69, 8.7))
plt.plot(x_pos, y_pos, marker='s', color='red', linewidth=1, markersize=3)
plt.xlabel('x / μm', fontsize=15)
plt.ylabel('y / μm', fontsize=15)
plt.title('Diagram 3.1: Brownian Motion Trajectory', fontsize=15)
plt.grid(True, alpha=0.3)
plt.axis('equal')
plt.savefig(f"Graph_Brown.pdf", format="pdf", bbox_inches='tight')
plt.show()


# %% Section 2: Calculate displacements and mean squared displacement
print("\n" + "="*50)
print("Section 2: Calculating displacements")
print("="*50)

# Calculate time intervals and position changes between consecutive points
dt = np.diff(time)  # Time differences
dx = np.diff(x_pos)  # x-position differences
dy = np.diff(y_pos)  # y-position differences

# Calculate squared displacement: r² = (Δx)² + (Δy)²
r_squared = dx**2 + dy**2

# Calculate statistics (excluding initial transient data points 0-23)
start_index = 23
r_squared_mean = np.mean(r_squared[start_index:])
r_squared_std = np.std(r_squared[start_index:]) / np.sqrt(len(r_squared[start_index:]))
dt_mean = np.mean(dt[start_index:])
dt_std = np.std(dt[start_index:]) / np.sqrt(len(dt[start_index:]))

print(f"Mean squared displacement: {r_squared_mean:.4f} ± {r_squared_std:.4f} μm²")
print(f"Mean time interval: {dt_mean:.4f} ± {dt_std:.6f} s")


# %% Section 3: Calculate Boltzmann constant and diffusion coefficient
print("\n" + "="*50)
print("Section 3: Physical constants calculation")
print("="*50)

# Calculate Boltzmann constant using Einstein relation
# k_B = (6πηa⟨r²⟩) / (4T⟨Δt⟩)
k_boltzmann = (6 * np.pi * PARTICLE_RADIUS * VISCOSITY * r_squared_mean * 1e-12) / \
              (4 * TEMPERATURE * dt_mean)

# Calculate uncertainty in Boltzmann constant
dk_boltzmann = k_boltzmann * np.sqrt(
    (TEMP_UNCERTAINTY/TEMPERATURE)**2 + 
    (VISCOSITY_UNCERTAINTY/VISCOSITY)**2 + 
    (RADIUS_UNCERTAINTY/PARTICLE_RADIUS)**2 + 
    (r_squared_std/r_squared_mean)**2 + 
    (dt_std/dt_mean)**2
)

# Calculate diffusion coefficient: D = k_B·T / (6πηa)
D = k_boltzmann * TEMPERATURE / (6 * np.pi * VISCOSITY * PARTICLE_RADIUS)

# Calculate uncertainty in diffusion coefficient
dD = D * np.sqrt(
    (dk_boltzmann/k_boltzmann)**2 + 
    (TEMP_UNCERTAINTY/TEMPERATURE)**2 + 
    (VISCOSITY_UNCERTAINTY/VISCOSITY)**2 + 
    (RADIUS_UNCERTAINTY/PARTICLE_RADIUS)**2
)

print(f"\nResults:")
print(f"Boltzmann constant: ({k_boltzmann:.4e} ± {dk_boltzmann:.4e}) J/K")
print(f"Literature value:    1.381e-23 J/K")
print(f"Diffusion coefficient: ({D:.4e} ± {dD:.4e}) m²/s")


# %% Section 4: Histogram of displacements
print("\n" + "="*50)
print("Section 4: Creating displacement histogram")
print("="*50)

# Import scipy.stats for the normal distribution
from scipy.stats import norm

# Combine x and y displacements for histogram
all_displacements = np.append(dx[start_index:-1], dy[start_index:-1])

# Calculate statistics for Gaussian fit
mu = np.mean(all_displacements)
sigma = np.std(all_displacements)

print(f"Displacement distribution: μ = {mu:.2f} μm, σ = {sigma:.2f} μm")

# Create histogram
plt.figure(figsize=(11.69, 8.7))
n, bins, patches = plt.hist(
    all_displacements, 
    bins='auto', 
    density=True, 
    color='r', 
    ec='black', 
    alpha=0.7
)

# Overlay Gaussian distribution using scipy.stats.norm
x_range = np.linspace(all_displacements.min(), all_displacements.max(), 100)
gaussian = norm.pdf(x_range, mu, sigma)
plt.plot(x_range, gaussian, 'b-', linewidth=2, label='Gaussian fit')

# Add statistics text
plt.text(
    0.7, 0.85, 
    f'μ = {mu:.2f} μm\nσ = {sigma:.2f} μm', 
    transform=plt.gca().transAxes,
    fontsize=20, 
    verticalalignment='top',
    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5)
)

plt.xlabel('Displacement / μm', fontsize=20)
plt.ylabel('Relative Frequency', fontsize=20)
plt.title('Diagram 3.2: Histogram of Displacements', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(color='black', linestyle='-.', alpha=0.3)
plt.legend(fontsize=15)
plt.savefig(f"Brown_histogram.pdf", format="pdf", bbox_inches='tight')
plt.show()

# %% Section 5: Cumulative squared displacement
print("\n" + "="*50)
print("Section 5: Cumulative displacement analysis")
print("="*50)

# Define linear function for fitting
def linear(x, slope, intercept):
    """Linear function y = slope*x + intercept"""
    return slope * x + intercept

# Calculate cumulative sum of squared displacements
r_cumulative = np.cumsum(r_squared)

# Note: r_cumulative has length len(time)-1 because r_squared comes from diff()
# We need to use time[1:] to match the length of r_cumulative
time_for_cumulative = time[1:]  # Skip first time point to match r_cumulative length

# Fit linear function to cumulative data (excluding initial transient)
popt, pcov = curve_fit(
    linear, 
    time_for_cumulative[start_index:], 
    r_cumulative[start_index:]
)
slope, intercept = popt
slope_err = np.sqrt(pcov[0, 0])

print(f"Linear fit: slope = {slope:.4f} ± {slope_err:.4f} μm²/s")
print(f"Expected slope from Einstein relation: {4*D*1e12:.4f} μm²/s")

# Plot cumulative displacement
plt.figure(figsize=(11.69, 8.7))
plt.plot(
    time_for_cumulative[start_index:], 
    r_cumulative[start_index:], 
    marker='.', 
    markersize=2, 
    color='red', 
    linewidth=0,
    label='Experimental data'
)
plt.plot(
    time_for_cumulative[start_index:], 
    linear(time_for_cumulative[start_index:], *popt), 
    'b-', 
    linewidth=2,
    label=f'Linear fit (slope={slope:.2f} μm²/s)'
)

plt.xlabel('Time / s', fontsize=15)
plt.ylabel('Cumulative Σr² / μm²', fontsize=15)
plt.title('Diagram 3.3: Cumulative Squared Displacement', fontsize=15)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3)
plt.savefig(f"Brown_cumulative.pdf", format="pdf", bbox_inches='tight')
plt.show()

print("\n" + "="*50)
print("Analysis complete!")
print("="*50)