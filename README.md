# Journal Bearing Hydrodynamics Solver

## Overview
This repository contains a Python solver that numerically evaluates the 2D Reynolds equation for a finite journal bearing. By utilizing the Finite Difference Method (FDM) and a Gauss-Seidel iterative technique, the script calculates the non-dimensional pressure distribution across the lubricating fluid film. 

Crucially, the numerical model applies the Reynolds (Half-Sommerfeld) boundary condition to truncate negative pressures, accurately simulating real-world fluid cavitation in the diverging clearance of the bearing.

### Key Features
* **Pressure Field Calculation:** Solves for non-dimensional pressure over a discrete M x N spatial grid.
* **Data Visualization:** Automatically plots the circumferential pressure profile at the bearing centerplane (z = 0.5).
* **Macroscopic Parameter Extraction:** Utilizes Simpson's 1/3rd Rule (via SciPy) to integrate the pressure field, outputting the total load-carrying capacity and attitude angle across various L/D and eccentricity ratios.

## Theory & Mathematical Derivation
For a complete discussion of the governing physical equations, the step-by-step finite difference discretization, and a formal analysis of the results, please see the [Project Report](docs/report.pdf) located in the `docs/` directory.

## Prerequisites
To run this code, you will need Python installed on your system along with the following standard scientific computing libraries:
```bash
pip install numpy matplotlib scipy
```

## How to Run
1. Clone this repository to your local machine.
2. Open your terminal and navigate to the project directory.
3. Execute the Python script by running: 
```bash
python main.py
```
*Note: Upon execution, a 2D graph will display the pressure distribution for L/D = 0.5 and e = 0.9. **You must close the plot window** for the script to continue and generate the performance tables in your terminal.*

## References
1. Lund, J. W., & Thomsen, K. K. (1982). *A Calculation Method and Data for the Dynamic Coefficients of Oil-Lubricated Journal Bearings*. Topics in Fluid Film Bearing and Rotor Bearing System Design and Optimization, ASME.
