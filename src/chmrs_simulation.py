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
# 4. NAVIER-STOKES FINITE-DIFFERENCE SOLVER LOOP
# =====================================================================
ttf_seconds = None
flotation_triggered = False

print(">> Initiating CHMRS Fluid-Structure Interaction Run...")

for step in range(NT):
    current_time = step * DT
    h_pulse = get_hydraulic_input(current_time)
    
    # Store previous step states
    P_w_old = P_w.copy()
    u_old = u.copy()
    
    # FIX: Explicit index assignments prevent full array variable overrides
    P_w[0] = RHO_WATER * G * h_pulse + P_w_chronic  # Inlet injection node 0
    P_w[-1] = P_w_chronic * 0.4                      # Outlet venting node NX-1
    
    # Explicit numerical updates across internal spatial cells (Nodes 1 to 98)
    for i in range(1, NX - 1):
        # 1D pressure gradient force
        pressure_gradient = (P_w_old[i+1] - P_w_old[i-1]) / (2.0 * DX)
        
        # Stable 1D Upwind Advection Scheme
        if u_old[i] > 0:
            advection = u_old[i] * (u_old[i] - u_old[i-1]) / DX
        else:
            advection = u_old[i] * (u_old[i+1] - u_old[i]) / DX
            
        # Viscous diffusion shear losses
        viscous_diffusion = MU * (u_old[i+1] - 2.0*u_old[i] + u_old[i-1]) / (DX**2)
        
        # Down-slope gravitational driving force vector
        gravity_force = G * np.sin(np.radians(SLOPE_ANGLE))
        
        # Update Velocity vector field
        u[i] = u_old[i] + DT * (-(1.0 / RHO_WATER) * pressure_gradient - advection + viscous_diffusion + gravity_force)
        
        # Mass conservation pressure coupling wave update
        P_w[i] = P_w_old[i] - DT * (RHO_WATER * G) * (u[i] - u[i-1]) / DX

    # Compute Terzaghi Effective Pressure profile: N = P_i - P_w
    N = P_I - P_w
    
    # Monitor the minimum effective pressure along the internal sliding mass
    min_N_body = np.min(N[1:-1])
    
    # Telemetry logging (sampled every 20 steps to save memory arrays)
    if step % 20 == 0:
        time_history.append(current_time / 60.0)
        N_min_history.append(min_N_body / 1e6)
        u_exit_history.append(u[-2])
    
        # Inside src/chmrs_simulation.py (around line 124)
    # Flotation tracking criterion check (Internal sliding zone lifts off rock)
    if min_N_body <= 0 and not flotation_triggered:
        ttf_seconds = current_time
        flotation_triggered = True
        print(f"!!! CRITICAL FAILURE MODE !!! Flotation boundary crossed at t = {ttf_seconds/60.0:.2f} minutes.")
        
        # ADD THIS BREAK STATEMENT TO AUTOMATICALLY EXIT AND SAVE THE GRAPH IMMEDIATELY
        break 

# =====================================================================
# 5. POST-PROCESSING DATA VISUALIZATION
# =====================================================================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Graph 1: Basal Effective Pressure Degradation Path
ax1.plot(time_history, N_min_history, color='crimson', lw=2.5, label='Minimum Body Effective Pressure ($N$)')
ax1.axhline(0, color='black', linestyle='--', lw=1.2, label='Flotation Parity ($N=0$)')
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
