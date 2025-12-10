import numpy as np
from parameters import *


def apply_boundary_conditions(lb, u_in=0.05):
    apply_periodic_x(lb)
    apply_bounce_back_walls(lb)
    apply_corner_bounce_back(lb)


def apply_periodic_x(lb):
    """Apply periodic boundary conditions in x-direction"""
    # This ensures what leaves on right enters on left, and vice versa

    # Method 1: Simple and effective
    for i in range(Q):
        # Copy right boundary (x=nx-1) from left interior (x=0)
        lb.f[i, -1, :] = lb.f[i, 0, :]

        # Copy left boundary (x=0) from right interior (x=nx-1)
        lb.f[i, 0, :] = lb.f[i, -1, :]


#  BOUNCE-BACK WALLS (TOP AND BOTTOM)
def apply_bounce_back_walls(lb):
    """
    applying standard bounce-back for a no-slip channel.
    """

    # Bottom wall (y = 0)
    y = 0
    lb.f[2, :, y] = lb.f[4, :, y]
    lb.f[5, :, y] = lb.f[7, :, y]
    lb.f[6, :, y] = lb.f[8, :, y]

    # Top wall (y = ny-1)
    y = lb.ny - 1
    lb.f[4, :, y] = lb.f[2, :, y]
    lb.f[7, :, y] = lb.f[5, :, y]
    lb.f[8, :, y] = lb.f[6, :, y]


# CORNERS
def apply_corner_bounce_back(lb):
    """
    Applying pure bounce-back at the four corners.
    """

    # bottom-left corner (0,0)
    lb.f[2, 0, 0] = lb.f[4, 0, 0]
    lb.f[5, 0, 0] = lb.f[7, 0, 0]
    lb.f[6, 0, 0] = lb.f[8, 0, 0]

    # top-left corner (0, ny-1)
    y = lb.ny - 1
    lb.f[4, 0, y] = lb.f[2, 0, y]
    lb.f[7, 0, y] = lb.f[5, 0, y]
    lb.f[8, 0, y] = lb.f[6, 0, y]

    # bottom-right corner (nx-1, 0)
    x = lb.nx - 1
    lb.f[2, x, 0] = lb.f[4, x, 0]
    lb.f[5, x, 0] = lb.f[7, x, 0]
    lb.f[6, x, 0] = lb.f[8, x, 0]

    # top-right corner (nx-1, ny-1)
    lb.f[4, x, y] = lb.f[2, x, y]
    lb.f[7, x, y] = lb.f[5, x, y]
    lb.f[8, x, y] = lb.f[6, x, y]
