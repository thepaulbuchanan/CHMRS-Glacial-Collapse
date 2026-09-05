"""
CHMRS Numerical Simulation Engine (Paper 3 Companion Code)
Author: Paul Buchanan
Date: September 5, 2026

Description:
Solves a 1D coupled Navier-Stokes fluid pressure wave under a Heaviside 
episodic rainfall shock pulse (e.g., 50mm over 2 hours) to determine the 
basal Effective Pressure (N) evolution and identify the Time-to-Failure (TTF) 
flotation threshold (N <= 0).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq  # ADD THIS LINE FOR SPECTRAL MATH


# =====================================================================
# 1. PHYSICAL CONSTANTS & CONTEXT PRESETS
# =====================================================================
RHO_ICE = 917.0      # Density of glacial ice (kg/m^3)
RHO_WATER = 1000.0   # Density of fresh water (kg/m^3)
G = 9.81             # Acceleration due to gravity (m/s^2)
MU = 1.002e-3        # Dynamic viscosity of water at 20C (Pa*s)

ICE_THICKNESS = 120.0  # Core thickness H (m)
SLOPE_ANGLE = 32.0     # Steep hanging glacier incline (degrees)
BED_LENGTH = 1500.0    # Horizontal reach of the sliding zone (m)

P_I = RHO_ICE * G * ICE_THICKNESS  # Static Ice Overburden Pressure (Pa)

# =====================================================================
# 2. NUMERICAL GRID AND TIME DISCRETIZATION
# =====================================================================
NX = 100               # Number of spatial grid points along the bedrock line
DX = BED_LENGTH / NX   # Spatial step size (m)
NT = 144000            # Total steps for a 2-hour window (144000 * 0.05s = 7200s)
DT = 0.05              # Temporal step size (s) - strictly satisfies CFL stability

x = np.linspace(0, BED_LENGTH, NX)
t = np.linspace(0, NT * DT, NT)

P_w = np.zeros(NX)     # Basal Water Pressure field (Pa)
u = np.zeros(NX)       # Subglacial fluid velocity vector (m/s)

# Initialize system to a stable, unfloated baseline pressure
P_w_chronic = 0.40 * P_I  
P_w[:] = np.linspace(P_w_chronic, P_w_chronic * 0.4, NX)

# =====================================================================
# 3. HEAVISIDE EPISODIC SHOCK PULSE FORCING
# =====================================================================
def get_hydraulic_input(current_time):
    """Implements the bounded Heaviside step function for rainfall pulse."""
    t_start = 0.0
    t_end = 7200.0
    if t_start <= current_time <= t_end:
        return 0.94 * ICE_THICKNESS  # Triggers flotation conditions
    else:
        return 0.0

time_history = []
N_min_history = []
u_exit_history = []

# =====================================================================
# 4. NAVIER-STOKES SOLVER WITH GEOMETRIC NOISE & STOCHASTIC LÉVY JUMPS
# =====================================================================
ttf_seconds = None
flotation_triggered = False

# Configuration parameters for stochastic boundary layer
SIGMA_GEOM = 0.08 * P_I  # Geometric noise intensity (8% of ice weight)
np.random.seed(42)       # Fix seed for strict reproducibility

# --- NEW: STOCHASTIC LÉVY JUMP ARCHITECTURE PARAMETERS ---
LAMBDA_JUMP = 0.005      # Probability of a structural choke/jump occurring per step (0.5% chance)
GAMMA_JUMP = 0.15 * P_I  # Scale amplitude of a structural pressure shock wave jump

prob_flotation_triggered = False
ttf_prob_seconds = None

# Allocate the grid array space cleanly across the 100 spatial nodes
P_w = np.linspace(P_w_chronic, P_w_chronic * 0.4, NX)

print(">> Initiating Coupled Stochastic Lévy Jump Boundary Run...")

for step in range(NT):
    current_time = step * DT
    h_pulse = get_hydraulic_input(current_time)
    
    # Store previous step states
    P_w_old = P_w.copy()
    u_old = u.copy()
    
    # Explicit boundary conditions to upstream/downstream terminal nodes
    # Inside src/chmrs_simulation.py -> Line 91 (Fix the variable overwrite)
    P_w[0] = RHO_WATER * G * h_pulse + P_w_chronic  # Target index node 0 exclusively

    P_w[-1] = P_w_chronic * 0.4                     # Outlet venting node NX-1
    
    # --- DETERMINISTIC NAVIER-STOKES SOLVER STEP ---
    for i in range(1, NX - 1):
        pressure_gradient = (P_w_old[i+1] - P_w_old[i-1]) / (2.0 * DX)
        
        if u_old[i] > 0:
            advection = u_old[i] * (u_old[i] - u_old[i-1]) / DX
        else:
            advection = u_old[i] * (u_old[i+1] - u_old[i]) / DX
            
        viscous_diffusion = MU * (u_old[i+1] - 2.0*u_old[i] + u_old[i-1]) / (DX**2)
        gravity_force = G * np.sin(np.radians(SLOPE_ANGLE))
        
        u[i] = u_old[i] + DT * (-(1.0 / RHO_WATER) * pressure_gradient - advection + viscous_diffusion + gravity_force)
        P_w[i] = P_w_old[i] - DT * (RHO_WATER * G) * (u[i] - u[i-1]) / DX

    # --- PROBABILISTIC GEOMETRIC NOISE & STOCHASTIC LÉVY JUMP LAYER ---
    # 1. Continuous Gaussian noise field
    geometric_noise = np.random.normal(0, SIGMA_GEOM, size=NX)
    
    # 2. Compound Poisson Jump Process (Lévy Jump Matrix)
    # Checks if a sudden subglacial structural choking event triggers at this millisecond step
    levy_jump = np.zeros(NX)
    if np.random.rand() < LAMBDA_JUMP:
        # A roof collapse occurs! Generate a massive, positive pressure spike spike along the bed
        jump_location = np.random.randint(1, NX - 1)
        levy_jump[jump_location] = np.random.normal(GAMMA_JUMP, 0.05 * GAMMA_JUMP)
        print(f"   [LÉVY SHOCK] Subglacial conduit collapse modeled at cell node {jump_location} at t = {current_time/60.0:.2f} mins.")
    
    # Compute the decoupled probabilistic pressure field containing discontinuous shocks
    P_w_prob = P_w + geometric_noise + levy_jump
    
    # Compute Terzaghi Effective Pressure profiles: N = Pi - Pw
    N_deterministic = P_I - P_w
    N_probabilistic = P_I - P_w_prob
    
    min_N_det_body = np.min(N_deterministic[1:-1])
    min_N_prob_body = np.min(N_probabilistic[1:-1])
    
    # Telemetry logging (sampled every 20 steps to save memory arrays)
    if step % 20 == 0:
        time_history.append(current_time / 60.0)
        N_min_history.append(min_N_det_body / 1e6)
        u_exit_history.append(u[-2])
    
    # Evaluate Noise-Induced Transition Threshold (Probabilistic Localised Flotation)
    if min_N_prob_body <= 0 and not prob_flotation_triggered:
        ttf_prob_seconds = current_time
        prob_flotation_triggered = True
        print(f"--> PROBABILISTIC ALARM: Localised noise-induced flotation at t = {ttf_prob_seconds/60.0:.2f} minutes.")
        
    # Evaluate Deterministic Mean Flotation Limit (Wholesale Flotation)
    if min_N_det_body <= 0 and not flotation_triggered:
        ttf_seconds = current_time
        flotation_triggered = True
        print(f"!!! DETERMINISTIC FAILURE: Wholesale flotation at t = {ttf_seconds/60.0:.2f} minutes.")
        
    # AUTOMATED EARLY EXIT: Once both boundaries are triggered, compute the delta and break
    if flotation_triggered and prob_flotation_triggered:
        delta_lead_time = (ttf_seconds - ttf_prob_seconds) / 60.0
        print(f">> SUCCESS: Stochastic Lévy-diffusion delta identifies localized failure {delta_lead_time:.2f} minutes early.")
        break

# =====================================================================
# 5. POST-PROCESSING DATA VISUALIZATION & SPECTRAL FFT MODULE
# =====================================================================
# Convert telemetry to clean numpy arrays
t_arr = np.array(time_history)        # Time in minutes
N_arr = np.array(N_min_history)       # Effective Pressure in MPa
u_arr = np.array(u_exit_history)      # Venting velocity in m/s

# Calculate the actual sampling rate for the FFT computation
# Sampling interval in seconds: 20 steps * DT (0.05s) = 1.0 second
dt_sample = 20 * DT  
N_samples = len(u_arr)

# Execute the Fast Fourier Transform over the velocity profile
u_detrended = u_arr - np.mean(u_arr)  # Strip DC offset to cleanly see frequencies
fft_values = fft(u_detrended)
frequencies = fftfreq(N_samples, d=dt_sample)

# Calculate Power Spectral Density (PSD) - positive frequencies only
positive_mask = frequencies > 0
freq_spectrum = frequencies[positive_mask]
psd_spectrum = np.abs(fft_values[positive_mask])**2

# --- COMPILE PEER-REVIEW GRID VISUALIZATION ---
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 11))

# Graph 1: Basal Effective Pressure Degradation Path
ax1.plot(t_arr, N_arr, color='crimson', lw=2.5, label='Min Body Effective Pressure ($N$)')
ax1.axhline(0, color='black', linestyle='--', lw=1.2, label='Flotation Threshold ($N=0$)')
if ttf_prob_seconds is not None:
    ax1.axvline(ttf_prob_seconds/60.0, color='forestgreen', linestyle=':', lw=2,
                label=f'Stochastic Alarm: {ttf_prob_seconds/60.0:.2f} mins')
if ttf_seconds is not None:
    ax1.axvline(ttf_seconds/60.0, color='darkorange', linestyle=':', lw=2,
                label=f'Wholesale Failure: {ttf_seconds/60.0:.2f} mins')
ax1.set_ylabel('Effective Pressure $N$ (MPa)', fontsize=11)
ax1.set_title('CHMRS Solver: Coupled Stochastic Boundary Layer & Spectral Output', fontsize=13, fontweight='bold')
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(loc='upper right')

# Graph 2: Forward-Axis Conduit Venting Velocity Pulse
ax2.plot(t_arr, u_arr, color='dodgerblue', lw=2.5, label='Terminal Venting Velocity ($u$)')
ax2.set_ylabel('Fluid Outflow $u$ (m/s)', fontsize=11)
ax2.grid(True, linestyle=':', alpha=0.6)
ax2.legend(loc='upper left')

# Graph 3: Fast Fourier Transform Power Spectral Density (The Spectral Blueprint)
ax3.plot(freq_spectrum, psd_spectrum, color='purple', lw=2.5, label='Velocity Power Spectrum (PSD)')
ax3.set_xlabel('Frequency (Hz)', fontsize=11)
ax3.set_ylabel('Power Density $(m/s)^2/Hz$', fontsize=11)
ax3.set_title('Spectral Signature: Critical Slowing Down Dominance Pattern', fontsize=11, fontweight='bold')
ax3.set_yscale('log')  # Log scale reveals the low-frequency resonance spike cleanly
ax3.grid(True, linestyle=':', alpha=0.6)
ax3.legend(loc='upper right')

plt.tight_layout()
plt.savefig('chmrs_simulation_output.png', dpi=300)
print(">> Spectral analytics plotted successfully. Exported to: chmrs_simulation_output.png")
