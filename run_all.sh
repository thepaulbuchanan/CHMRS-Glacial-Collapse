#!/bin/bash
# CHMRS Automated Execution and Validation Shell Script
# Author: Paul Buchanan
# Date: September 5, 2026

echo "========================================================="
# 1. Fire up the programmatic unit tests first as a gatekeeper
echo ">> Step 1: Running Navier-Stokes CFL Stability Verification..."
python3 -m unittest tests/test_simulation.py

# Capture the exit code of the test engine
if [ $? -eq 0 ]; then
    echo ">> [PASS] Grid parameters verified. Proceeding to execution."
    echo "========================================================="
    
    # 2. Run the main solver to compile your data visualization chart
    echo ">> Step 2: Launching CHMRS Non-Linear Simulation Solver..."
    python3 src/chmrs_simulation.py
else
    echo ">> [FAIL] Code stability checks failed. Halting runtime loop."
    exit 1
fi
