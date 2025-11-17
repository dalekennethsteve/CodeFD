import numpy as np
from parameters import *


def apply_boundary_conditions(lb):
    """Apply simplified boundary conditions"""

    # Inlet (left boundary) - constant velocity profile
    apply_inlet_velocity(lb)

    # Outlet (right boundary) - constant pressure
    apply_outlet_pressure(lb)

    # Walls (top and bottom)
    apply_wall_bounce_back(lb)


def apply_inlet_velocity(lb):
    """Apply constant velocity at inlet"""
    x = 0  # left boundary
    u_in = 0.05  # inlet velocity

    # Parabolic velocity profile for Poiseuille flow
    for y in range(lb.ny):
        if y == 0 or y == lb.ny - 1:  # walls
            u_x = 0.0
        else:
            # Parabolic profile
            u_x = u_in * 1.5 * (1 - ((y - lb.ny / 2) / (lb.ny / 2)) ** 2)

        # Simple velocity boundary condition
        lb.u[0, x, y] = u_x
        lb.u[1, x, y] = 0.0
        lb.rho[x, y] = rho0

        # Recalculate equilibrium at boundary
        for i in range(Q):
            cu = c[i, 0] * lb.u[0, x, y] + c[i, 1] * lb.u[1, x, y]
            u_sqr = lb.u[0, x, y] ** 2 + lb.u[1, x, y] ** 2
            lb.f[i, x, y] = w[i] * lb.rho[x, y] * (1 + 3 * cu + 4.5 * cu ** 2 - 1.5 * u_sqr)


def apply_outlet_pressure(lb):
    """Apply constant pressure at outlet"""
    x = lb.nx - 1  # right boundary

    for y in range(lb.ny):
        lb.rho[x, y] = rho0
        # Extrapolate velocity from interior
        if x > 0:
            lb.u[0, x, y] = lb.u[0, x - 1, y]
            lb.u[1, x, y] = lb.u[1, x - 1, y]

        # Recalculate equilibrium at boundary
        for i in range(Q):
            cu = c[i, 0] * lb.u[0, x, y] + c[i, 1] * lb.u[1, x, y]
            u_sqr = lb.u[0, x, y] ** 2 + lb.u[1, x, y] ** 2
            lb.f[i, x, y] = w[i] * lb.rho[x, y] * (1 + 3 * cu + 4.5 * cu ** 2 - 1.5 * u_sqr)


def apply_wall_bounce_back(lb):
    """Apply bounce-back boundary conditions at top and bottom walls"""
    # Bottom wall (y = 0)
    for x in range(lb.nx):
        # Store populations that hit the wall
        lb.f[2, x, 0] = lb.f[4, x, 0]  # bounce back from bottom
        lb.f[5, x, 0] = lb.f[7, x, 0]
        lb.f[6, x, 0] = lb.f[8, x, 0]

        # Set wall velocity to zero
        lb.u[0, x, 0] = 0.0
        lb.u[1, x, 0] = 0.0

    # Top wall (y = lb.ny-1)
    for x in range(lb.nx):
        lb.f[4, x, -1] = lb.f[2, x, -1]  # bounce back from top
        lb.f[7, x, -1] = lb.f[5, x, -1]
        lb.f[8, x, -1] = lb.f[6, x, -1]

        # Set wall velocity to zero
        lb.u[0, x, -1] = 0.0
        lb.u[1, x, -1] = 0.0