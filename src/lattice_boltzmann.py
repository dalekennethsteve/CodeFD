import numpy as np
from parameters import *
from boundary_conditions import apply_boundary_conditions

class LatticeBoltzmann:
    def __init__(self):
        self.nx = nx
        self.ny = ny
        self.f = np.zeros((Q, nx, ny))
        self.f_eq = np.zeros((Q, nx, ny))
        self.rho = np.ones((nx, ny)) * rho0
        self.u = np.zeros((2, nx, ny))
        self.initialize()

    def initialize(self):
        """Initialize with equilibrium"""
        self.rho[:, :] = rho0
        self.u[0, :, :] = 0.0
        self.u[1, :, :] = 0.0
        self.calculate_equilibrium()
        self.f[:, :, :] = self.f_eq[:, :, :]

    def calculate_equilibrium(self):
        """Calculate equilibrium distribution"""
        u_sqr = self.u[0]**2 + self.u[1]**2
        for i in range(Q):
            cu = c[i,0]*self.u[0] + c[i,1]*self.u[1]
            self.f_eq[i] = w[i] * self.rho * (1 + 3*cu + 4.5*cu**2 - 1.5*u_sqr)

    def apply_body_force_guo(self):
        Fx, Fy = body_force
        for i in range(Q):
            ciux = c[i, 0] - self.u[0]
            ciuy = c[i, 1] - self.u[1]
            cu = c[i, 0] * self.u[0] + c[i, 1] * self.u[1]
            Fi = ((ciux + 3 * cu * c[i, 0]) * Fx +
                  (ciuy + 3 * cu * c[i, 1]) * Fy)
            self.f[i] += w[i] * (1 - 0.5 * omega) * 3.0 * Fi

    def stream(self):
        """Simple streaming with numpy rolls"""
        f_temp = self.f.copy()
        for i in range(Q):
            # Use numpy roll for efficiency
            self.f[i] = np.roll(f_temp[i], shift=(c[i,0], c[i,1]), axis=(0,1))

    def collide(self):
        """Collision step - CORRECT ORDER"""
        self.calculate_macroscopic()  # Get ρ,u from current f
        self.calculate_equilibrium()  # Calculate f_eq
        self.apply_body_force_guo()   # Apply force BEFORE collision
        self.f -= omega * (self.f - self.f_eq)  # Collision

    def calculate_macroscopic(self):
        self.rho = np.sum(self.f, axis=0)
        mom_x = np.sum(self.f * c[:, 0].reshape(Q, 1, 1), axis=0)
        mom_y = np.sum(self.f * c[:, 1].reshape(Q, 1, 1), axis=0)

        # Guo velocity correction
        self.u[0] = (mom_x + 0.5 * body_force[0]) / self.rho
        self.u[1] = (mom_y + 0.5 * body_force[1]) / self.rho



    def step(self):
        """ LBM sequence"""
        self.stream()
        apply_boundary_conditions(self)
        self.collide()