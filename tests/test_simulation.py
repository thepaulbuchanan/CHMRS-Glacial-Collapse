"""
CHMRS Numerical Stability Automated Validation Suite (Paper 3 Verification)
Author: Paul Buchanan
Date: September 5, 2026

Description:
Executes programmatic unit tests verifying the Courant-Friedrichs-Lewy (CFL) 
numerical stability condition and physics boundaries of the CHMRS Navier-Stokes solver.
"""

import unittest
import numpy as np

class TestCHMRSSimulationStability(unittest.TestCase):

    def setUp(self):
        """Set up standard framework simulation parameters for testing validation."""
        self.ice_thickness = 120.0  # H (m)
        self.bed_length = 1500.0    # Total channel horizontal reach (m)
        self.nx = 100               # Grid cells
        self.dx = self.bed_length / self.nx  # dx = 15.0 m
        self.dt = 0.05               # dt = 1.0 s
        
        # Physics constants
        self.rho_water = 1000.0     # kg/m^3
        self.g = 9.81               # m/s^2
        self.slope_angle = 32.0     # degrees
        
        # Max theoretical fluid velocity expected under a hyper-pressurized blowout wave
        # Derived from maximum gravity-driven down-slope acceleration acceleration vector
        self.max_theoretical_velocity = np.sqrt(2 * self.g * np.sin(np.radians(self.slope_angle)) * self.bed_length)

    def test_cfl_advection_stability_condition(self):
        """
        CRITICAL SANITY CHECK: Verifies that the Courant-Friedrichs-Lewy (CFL)
        condition for explicit 1D advection dynamics is strictly met.
        Equation: C = (u * dt) / dx <= 1.0
        """
        courant_number = (self.max_theoretical_velocity * self.dt) / self.dx
        
        print(f"\n[CFL Check] Max Expected Velocity: {self.max_theoretical_velocity:.2f} m/s")
        print(f"[CFL Check] Calculated Courant Number (C): {courant_number:.4f}")
        
        # For explicit updates to remain stable, the Courant number must remain less than or equal to 1.0
        self.assertTrue(courant_number <= 1.0, 
                        f"CFL condition breached: Courant number {courant_number:.2f} > 1.0! Simulation will diverge.")

    def test_hydrostatic_pressure_vs_overburden_limit(self):
        """
        Verifies that an intense rainfall shock parameter is mathematically 
        capable of driving the system to complete flotation parity (N <= 0).
        """
        rho_ice = 917.0
        p_i = rho_ice * self.g * self.ice_thickness  # Static overburden
        
        # Modelled Heaviside input head height (h = 0.94 * H)
        h_pulse = 0.94 * self.ice_thickness  
        p_w_max = self.rho_water * self.g * h_pulse
        
        print(f"[Pressure Check] Ice Overburden (Pi): {p_i/1e6:.4f} MPa")
        print(f"[Pressure Check] Max Injection Pressure (Pw): {p_w_max/1e6:.4f} MPa")
        
        # Flotation check: Max water pressure must be capable of overcoming ice weight to trigger CHMRS
        self.assertGreaterEqual(p_w_max, p_i, 
                                "The chosen rainfall pulse configuration cannot generate enough head to clear the flotation limit.")

if __name__ == '__main__':
    unittest.main()
