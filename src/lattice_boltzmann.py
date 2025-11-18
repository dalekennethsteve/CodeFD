import numpy as np
from parameters import *
from src.boundary_conditions import apply_boundary_conditions


class LatticeBoltzmann:
    def __init__(self):
        self.nx = nx
        self.ny = ny

        # Distribution functions
        self.f = np.zeros((Q, nx, ny))
        self.f_eq = np.zeros((Q, nx, ny))

        # Macroscopic variables
        self.rho = np.ones((nx, ny)) * rho0
        self.u = np.zeros((2, nx, ny))

        self.initialize()

    def initialize(self):
        """Initialize distribution functions with equilibrium values"""
        self.rho[:, :] = rho0
        self.u[0, :, :] = 0.0
        self.u[1, :, :] = 0.0
        self.calculate_equilibrium()
        self.f[:, :, :] = self.f_eq[:, :, :]

    def calculate_equilibrium(self):
        """Calculate equilibrium distribution function"""
        u_sqr = self.u[0] ** 2 + self.u[1] ** 2

        for i in range(Q):
            cu = c[i, 0] * self.u[0] + c[i, 1] * self.u[1]
            self.f_eq[i] = w[i] * self.rho * (1 + 3 * cu + 4.5 * cu ** 2 - 1.5 * u_sqr)

    def apply_body_force_guo(self):
        """Guo et al. (2002) body force scheme"""
        for i in range(Q):
            cu = c[i, 0] * self.u[0] + c[i, 1] * self.u[1]

            # Guo force term
            force_term = w[i] * (1 - 0.5 * omega) * (3.0 * (c[i, 0] - self.u[0]) + 9.0 * cu * c[i, 0]) * body_force[0] + w[i] * (1 - 0.5 * omega) * (
                                 3.0 * (c[i, 1] - self.u[1]) + 9.0 * cu * c[i, 1]
                         ) * body_force[1]

            self.f[i] += force_term

    def stream(self):
        """Streaming step with periodic boundaries"""
        f_temp = self.f.copy()

        for i in range(Q):
            for x in range(self.nx):
                for y in range(self.ny):
                    # Apply periodic boundary conditions in x-direction
                    x_new = (x + c[i, 0]) % self.nx
                    # Apply bounce-back in y-direction (handled separately)
                    y_new = y + c[i, 1]

                    # Only stream if not hitting top/bottom walls
                    if 0 <= y_new < self.ny:
                        self.f[i, x_new, y_new] = f_temp[i, x, y]

    def collide(self):
        """Collision step with body force"""
        self.calculate_macroscopic()
        self.calculate_equilibrium()
        self.f = self.f - omega * (self.f - self.f_eq)
        self.apply_body_force_guo()  # Apply force AFTER collision

    def calculate_macroscopic(self):
        """Calculate macroscopic variables from distribution functions"""
        self.rho = np.sum(self.f, axis=0)
        self.u[0] = np.sum(self.f * c[:, 0].reshape(Q, 1, 1), axis=0) / self.rho
        self.u[1] = np.sum(self.f * c[:, 1].reshape(Q, 1, 1), axis=0) / self.rho

    def step(self):
        """Perform one LBM time step"""
        self.stream()
        apply_boundary_conditions(self)
        self.collide()