import numpy as np
from parameters import *


def apply_boundary_conditions(lb):
    """Apply periodic boundary conditions in x-direction, bounce-back in y-direction"""

    # Apply periodic boundary conditions in x-direction
    apply_periodic_x(lb)

    # Apply bounce-back at top and bottom walls
    apply_zou_he_walls(lb)



def apply_periodic_x(lb):
    """Apply periodic boundary conditions in x-direction"""
    # Copy distribution functions from right to left boundary and vice versa
    for i in range(Q):
        # Left boundary (x=0) gets populations from right boundary (x=nx-1)
        lb.f[i, 0, :] = lb.f[i, lb.nx - 1, :]

        # Right boundary (x=nx-1) gets populations from left boundary (x=0)
        lb.f[i, lb.nx - 1, :] = lb.f[i, 0, :]


def apply_zou_he_walls(lb):
    """Apply Zou-He boundary conditions at top and bottom walls"""

    # Bottom wall (y = 0) - no-slip
    for x in range(lb.nx):
        # Known distributions after streaming
        f0 = lb.f[0, x, 0]
        f1 = lb.f[1, x, 0]
        f3 = lb.f[3, x, 0]
        f4 = lb.f[4, x, 0]  # coming from wall
        f7 = lb.f[7, x, 0]
        f8 = lb.f[8, x, 0]

        # Wall velocity (no-slip)
        ux_wall = 0.0
        uy_wall = 0.0

        # Calculate density
        rho = (f0 + f1 + f3 + 2 * (f4 + f7 + f8)) / (1 - uy_wall)

        # Zou-He equations for unknown distributions
        lb.f[2, x, 0] = f4 + (2 / 3) * rho * uy_wall  # f2
        lb.f[5, x, 0] = f7 - 0.5 * (f1 - f3) + 0.5 * rho * ux_wall + (1 / 6) * rho * uy_wall  # f5
        lb.f[6, x, 0] = f8 + 0.5 * (f1 - f3) - 0.5 * rho * ux_wall + (1 / 6) * rho * uy_wall  # f6

        # Update macroscopic variables at boundary
        lb.rho[x, 0] = rho
        lb.u[0, x, 0] = ux_wall
        lb.u[1, x, 0] = uy_wall

    # Top wall (y = ny-1) - similar implementation
    for x in range(lb.nx):
        f0 = lb.f[0, x, -1]
        f1 = lb.f[1, x, -1]
        f2 = lb.f[2, x, -1]  # coming from wall
        f3 = lb.f[3, x, -1]
        f5 = lb.f[5, x, -1]
        f6 = lb.f[6, x, -1]

        ux_wall = 0.0
        uy_wall = 0.0

        rho = (f0 + f1 + f3 + 2 * (f2 + f5 + f6)) / (1 + uy_wall)  # note sign change

        lb.f[4, x, -1] = f2 - (2 / 3) * rho * uy_wall  # f4
        lb.f[7, x, -1] = f5 + 0.5 * (f1 - f3) - 0.5 * rho * ux_wall - (1 / 6) * rho * uy_wall  # f7
        lb.f[8, x, -1] = f6 - 0.5 * (f1 - f3) + 0.5 * rho * ux_wall - (1 / 6) * rho * uy_wall  # f8

        lb.rho[x, -1] = rho
        lb.u[0, x, -1] = ux_wall
        lb.u[1, x, -1] = uy_wall