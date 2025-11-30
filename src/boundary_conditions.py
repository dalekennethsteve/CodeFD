import numpy as np
from parameters import *


def apply_boundary_conditions(lb, u_in=0.05):
    apply_periodic_outlet(lb)
    apply_bounce_back_walls(lb)
    apply_zou_he_velocity_inlet(lb, u_in)
    apply_corner_bounce_back(lb)



#  ZOU–HE VELOCITY INLET (LEFT SIDE)

def apply_zou_he_velocity_inlet(lb, u_in):
    """
    applying Zou–He velocity inlet at left boundary (x = 0).
    """

    x = 0
    ux = u_in
    uy = 0.0

    for y in range(1, lb.ny - 1):
        f0 = lb.f[0, x, y]
        f2 = lb.f[2, x, y]
        f3 = lb.f[3, x, y]
        f4 = lb.f[4, x, y]
        f6 = lb.f[6, x, y]
        f7 = lb.f[7, x, y]

        # Density with Guo forcing correction
        rho = (f0 + f2 + f4 + 2.0*(f3 + f6 + f7) +
               (2.0 / 3.0) * (1.0 / omega) * body_force[0]) / (1.0 - ux)

        # Incoming unknowns: f1, f5, f8
        lb.f[1, x, y] = f3 + (2.0/3.0) * rho * ux
        lb.f[5, x, y] = f7 + 0.5*(f4 - f2) + (1.0/6.0)*rho*ux
        lb.f[8, x, y] = f6 + 0.5*(f2 - f4) + (1.0/6.0)*rho*ux

        lb.rho[x, y] = rho



#   PERIODIC OUTLET
def apply_periodic_outlet(lb):
    """
    Proper periodic boundary at x = nx - 1 (right outlet).
    """

    # outlet copies from node before inlet (nx-2),
    # inlet copies from node before outlet (1)
    lb.f[:, -1, :] = lb.f[:, 1, :]
    lb.f[:,  0, :] = lb.f[:, -2, :]



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
