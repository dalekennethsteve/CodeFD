CodeFD - Lattice Boltzmann Method for Fluid Flow Simulations

A custom Python implementation of the Lattice Boltzmann Method (LBM) for simulating 2D laminar fluid flows with various obstacle geometries.

The goal of this project is to explore the Lattice Boltzmann Method as a numerical tool for simulating laminar fluid flow using various boundaries. This is achieved by combining a theoretical study of the method with a practical numerical implementation in Python. Several benchmark flow configurations are investigated at different reynolds numbers, including Poiseuille flow and two-dimensional flow past a circular cylinder, rotating cylinder, assembly of cylinders and squares, in order to assess the accuracy, stability, and physical fidelity of the LBM in capturing key flow phenomena. 


Requirements:
- Python 3.7+
- NumPy
- Matplotlib

File Structure:

CodeFD/
├── parameters.py                 # Physical and numerical parameters
├── lattice_boltzmann.py          # Core LBM solver class
├── boundary_conditions.py        # Boundary condition handling
├── visualization.py              # Plotting utilities
├── utils.py                      # Convergence checking
├── main.py                       # Single cylinder simulation
├── main_multipleCylinders.py     # Multiple cylinders simulation
├── main_multipleSquares.py       # Multiple squares simulation
├── main_rotating.py              # Rotating cylinder simulation
├── poiseuille.py                 # Poiseuille flow validation
├── rotating_cylinder.py          # Rotating cylinder solver
├── visualization_rotating.py     # Rotating cylinder plots


Simulation Parameters:

Domain size : 300 × 62 lattice nodes 
Relaxation time (τ) : 0.8 
Kinematic viscosity (ν) : 0.1 
Convergence tolerance : 1×10⁻⁶  
Max time steps : 10,000 

Authors: 
Kenneth Steve Diyya (127799)
Sai Pranay Unnam (126722)

Supervisors: Prof. Dr. rer. nat. Björn Rüffer, Dr. Natalia Gorban, M.Sc. Tighana Wenge Basele

Institution: Bauhaus-Universität Weimar

License: 
This project is for academic purposes. Contact the authors for more information.
