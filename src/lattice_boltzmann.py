import numpy as np
from parameters import *


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
        # Initialize with small velocity to start flow
        self.u[0, :, :] = 0.01  # small initial x-velocity
        self.calculate_equilibrium()
        self.f[:, :, :] = self.f_eq[:, :, :]

    def calculate_equilibrium(self):
        """Calculate equilibrium distribution function"""
        u_sqr = self.u[0] ** 2 + self.u[1] ** 2

        for i in range(Q):
            cu = c[i, 0] * self.u[0] + c[i, 1] * self.u[1]
            self.f_eq[i] = w[i] * self.rho * (1 + 3 * cu + 4.5 * cu ** 2 - 1.5 * u_sqr)

    def stream(self):
        """Streaming step - fixed implementation"""
        f_temp = self.f.copy()

        for i in range(Q):
            # Apply periodic boundary conditions manually
            for x in range(self.nx):
                for y in range(self.ny):
                    x_new = (x + c[i, 0]) % self.nx
                    y_new = (y + c[i, 1]) % self.ny
                    self.f[i, x_new, y_new] = f_temp[i, x, y]

    def collide(self):
        """Collision step"""
        self.calculate_macroscopic()
        self.calculate_equilibrium()
        self.f = self.f - omega * (self.f - self.f_eq)

    def calculate_macroscopic(self):
        """Calculate macroscopic variables from distribution functions"""
        self.rho = np.sum(self.f, axis=0)
        self.u[0] = np.sum(self.f * c[:, 0].reshape(Q, 1, 1), axis=0) / self.rho
        self.u[1] = np.sum(self.f * c[:, 1].reshape(Q, 1, 1), axis=0) / self.rho

    def step(self):
        """Perform one LBM time step"""
        self.stream()
        self.collide()