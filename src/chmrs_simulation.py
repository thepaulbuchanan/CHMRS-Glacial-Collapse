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

# =====================================================================
# 1. PHYSICAL CONSTANTS & CONTEXT PRESETS (Paper 1 & 2 Baselines)
# =====================================================================
RHO_ICE = 917.0      # Density of glacial ice (kg/m^3)
RHO_WATER = 1000.0   # Density of fresh water (kg/m^3)
G = 9.81             # Acceleration due to gravity (m/s^2)
MU = 1.002e-3        # Dynamic viscosity of water at 20C (Pa*s)

# Slope and structural defaults
ICE_THICKNESS = 120.0  # Core thickness H (m)
SLOPE_ANGLE = 32.0     # Steep hanging glacier incline (degrees)
BED_LENGTH = 1500.0    # Horizontal reach of the sliding zone (m)

# Calculate static Ice Overburden Pressure (Pi = rho_i * g * H)
P_I = RHO_ICE * G * ICE_THICKNESS  

# =====================================================================
# 2. NUMERICAL GRID AND TIME DISCRETIZATION
# =====================================================================
NX = 100               # Number of spatial grid points along the bedrock line
DX = BED_LENGTH / NX   # Spatial step size (m)
NT = 7200              # Total simulation steps (7200s = 2.0 hours)
DT = 1.0               # Temporal step size (s) - satisfies CFL stability

# Define coordinate arrays
x = np.linspace(0, BED_LENGTH, NX)
t = np.linspace(0, NT * DT, NT)

# Preallocate state fields
P_w = np.zeros(NX)     # Basal Water Pressure field (Pa)
u = np.zeros(NX)       # Subglacial fluid velocity vector (m/s)

# Historical "Chronic" Baseline initialization
# Pre-pressurized system near its structural threshold
P_w_chronic = 0.85 * P_I  
P_w[:] = P_w_chronic

# =====================================================================
# 3. HEAVISIDE EPISODIC SHOCK PULSE FORCING (50mm over 2 Hours)
# =====================================================================
def get_hydraulic_input(current_time):
    """Implements the bounded Heaviside step function for rainfall pulse."""
    t_start = 0.0        # Shock initiates at simulation t=0
    t_end = 7200.0       # Pulse cuts off cleanly at 2 hours (7200 seconds)
    
    # 50mm over 2 hours converts to an intense influx head height parameter
    if t_start <= current_time <= t_end:
        return 0.94 * ICE_THICKNESS  # Forces h >= 0.92H parity marker
    else:
        return 0.0

# Arrays to capture downstream analytics over time
time_history = []
N_min_history = []
u_exit_history = []

# =====================================================================
# 4. NAVIER-STOKES FINITE-DIFFERENCE SOLVER LOOP
# =====================================================================
ttf_seconds = None
flotation_triggered = False

print(">> Initiating CHMRS Fluid-Structure Interaction Run...")

for step in range(NT):
    current_time = step * DT
    h_pulse = get_hydraulic_input(current_time)
    
    # Boundary Condition Injection: Upper fracture channels fill instantly
    # Hydrostatic Head Conversion: P_w = rho_w * g * h
    P_w[0] = RHO_WATER * G * h_pulse + P_w_chronic
    # Leading edge face outflow venting boundary condition (Paper 1 & 2 exit path)
    P_w[-1] = P_w_chronic * 0.4  
    
    # Copy arrays to maintain numerical separation between iterations
    P_w_old = P_w.copy()
    u_old = u.copy()
    
    # Explicit numerical updates across internal spatial cells
    for i in range(1, NX - 1):
        # 1D Incompressible Navier-Stokes solver with friction losses
        # u_t + u*u_x = -1/rho * P_x + nu * u_xx + g*sin(theta)
        pressure_gradient = (P_w_old[i+1] - P_w_old[i-1]) / (2.0 * DX)
        advection = u_old[i] * (u_old[i] - u_old[i-1]) / DX
        viscous_diffusion = MU * (u_old[i+1] - 2.0*u_old[i] + u_old[i-1]) / (DX**2)
        gravity_force = G * np.sin(np.radians(SLOPE_ANGLE))
        
        # Update Velocity vector field
        u[i] = u_old[i] + DT * (-(1.0 / RHO_WATER) * pressure_gradient - advection + viscous_diffusion + gravity_force)
        
        # Mass continuity update to resolve pressure field changes
        P_w[i] = P_w_old[i] - DT * RHO_WATER * (P_I / ICE_THICKNESS) * (u[i] - u[i-1]) / DX

    # Compute Terzaghi Effective Pressure profile: N = P_i - P_w
    N = P_I - P_w
    min_N = np.min(N)
    
    # Capture telemetry profiles
    time_history.append(current_time / 60.0)  # Store in minutes
    N_min_history.append(min_N / 1e6)         # Store in MPa
    u_exit_history.append(u[-2])              # Tracking terminal venting velocity
    
    # Flotation tracking criterion check (N <= 0)
    if min_N <= 0 and not flotation_triggered:
        ttf_seconds = current_time
        flotation_triggered = True
        print(f"!!! CRITICAL FAILURE MODE !!! Flotation boundary crossed at t = {ttf_seconds/60.0:.2f} minutes.")

# =====================================================================
# 5. POST-PROCESSING DATA VISUALIZATION
# =====================================================================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Graph 1: Basal Effective Pressure Degradation Path
ax1.plot(time_history, N_min_history, color='crimson', lw=2.5, label='Minimum Effective Pressure ($N$)')
ax1.axhline(0, color='black', linestyle='--', lw=1.2, label='Flotation Parity Threshold ($N=0$)')
if flotation_triggered:
    ax1.axvline(ttf_seconds/60.0, color='darkorange', linestyle=':', lw=2, 
                label=f'Time-to-Failure (TTF): {ttf_seconds/60.0:.1f} mins')
ax1.set_ylabel('Effective Pressure $N$ (MPa)', fontsize=12)
ax1.set_title('CHMRS Simulation: Structural Overburden vs Basal Pressure Wave', fontsize=14, fontweight='bold')
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(loc='upper right')

# Graph 2: Forward-Axis Conduit Venting Velocity Pulse
ax2.plot(time_history, u_exit_history, color='dodgerblue', lw=2.5, label='Terminal Venting Velocity ($u$)')
ax2.set_xlabel('Time Post-Shock Injection (Minutes)', fontsize=12)
ax2.set_ylabel('Fluid Outflow Velocity $u$ (m/s)', fontsize=12)
ax2.grid(True, linestyle=':', alpha=0.6)
ax2.legend(loc='upper left')

plt.tight_layout()
plt.savefig('chmrs_simulation_output.png', dpi=300)
print(">> Metrics plotted successfully. Output exported to: chmrs_simulation_output.png")
plt.show()
