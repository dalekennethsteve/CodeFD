import numpy as np
# Lattice parameters
D = 2  # dimensions
Q = 9  # D2Q9 lattice

# Lattice velocities
c = np.array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],[1, 1], [-1, 1], [-1, -1], [1, -1]])

# Lattice weights
w = np.array([
    4/9,  # stationary
    1/9, 1/9, 1/9, 1/9,  # cardinal directions
    1/36, 1/36, 1/36, 1/36  # diagonal directions
 ])

# Reverse directions for bounce-back
opposite = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])

# Simulation parameters
nx = 250    # length
ny = 50   # breadth
tau = 0.8   # relaxation time
rho0 = 1.0  # initial density
u_in = 0.05
# Body force parameter
body_force = np.array([5e-5, 0.0])

# Derived parameters
nu = (tau - 0.5) / 3.0  # kinematic viscosity
omega = 1.0 / tau       # relaxation time

# Calculate Reynolds number for information
H = ny - 2  # channel height (excluding walls)
Re = (u_in * H) / nu

print(f"Reynolds number: Re = {Re:.2f}")