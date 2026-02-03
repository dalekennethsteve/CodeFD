# parameters.py - Add cylinder parameters
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
nx = 300 # length
ny = 62   # breadth
tau = 0.8   # relaxation time
rho0 = 1.0  # initial density
# Body force parameter
body_force = np.array([0.00005, 0.0])

# Derived parameters
nu = 0.1     # kinematic viscosity (tau - 0.5) / 3.0
omega = 1.0 / tau       # relaxation time

H = ny - 2
F = body_force[0] / rho0

# Predicted maximum velocity:
u_max_predicted = (F * H**2) / (8 * nu)

# Reynolds number
H = ny - 2  # channel height
Re = (u_max_predicted * H) / nu

# Cylinder parameters
cylinder_radius = 6
cylinder_x = nx // 4
cylinder_y = ny // 2

cs = 1/np.sqrt(3)

'''
PARAMETERS FOR SIMULATION OF FLOW PAST MULTIPLE CYLINDERS AND SQUARES

# Square cylinder parameters
cylinder_width = 8
cylinder_type = 'square' 


# Multiple squares (x_center, y_center)
# Add as many squares as you like
cylinders = [
    (nx // 4,     ny // 2),        # Square 1
    (nx // 2,     ny // 2 + 10),   # Square 2
    (3 * nx // 4, ny // 2 - 10),   # Square 3
]
'''
