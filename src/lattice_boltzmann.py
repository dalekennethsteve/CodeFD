# lattice_boltzmann.py - COMPLETELY CORRECTED
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

        # Create obstacle mask for cylinder
        self.obstacle = np.zeros((nx, ny), dtype=bool)
        self.create_cylinder_mask()

        self.initialize()

    def create_cylinder_mask(self):
        """Create obstacle mask for multiple cylinders"""
        for cx, cy in cylinders:
            for x in range(self.nx):
                for y in range(self.ny):
                    dx = x - cx
                    dy = y - cy
                    if dx * dx + dy * dy <= cylinder_radius ** 2:
                        self.obstacle[x, y] = True

    def initialize(self):
        """Initialize all fields"""
        self.rho[:, :] = rho0
        self.u[0, :, :] = 0.0
        self.u[1, :, :] = 0.0
        self.calculate_equilibrium()
        self.f[:, :, :] = self.f_eq[:, :, :]

    def calculate_equilibrium(self):
        """Calculate equilibrium distribution"""
        u_sqr = self.u[0] ** 2 + self.u[1] ** 2
        for i in range(Q):
            cu = c[i, 0] * self.u[0] + c[i, 1] * self.u[1]
            self.f_eq[i] = w[i] * self.rho * (1 + 3 * cu + 4.5 * cu ** 2 - 1.5 * u_sqr)

    def apply_body_force(self):
        """Apply body force using Guo scheme - CORRECTED"""
        Fx, Fy = body_force

        # Calculate velocity magnitude for force adjustment
        u_mag = np.sqrt(self.u[0] ** 2 + self.u[1] ** 2)

        for i in range(Q):
            cu = c[i, 0] * self.u[0] + c[i, 1] * self.u[1]

            # CORRECT Guo force term (simplified form)
            force_term = w[i] * (1 - 0.5 * omega) * (
                    3.0 * (c[i, 0] - self.u[0]) + 9.0 * cu * c[i, 0]
            ) * Fx

            if Fy != 0:  # Only add if Fy is non-zero
                force_term += w[i] * (1 - 0.5 * omega) * (
                        3.0 * (c[i, 1] - self.u[1]) + 9.0 * cu * c[i, 1]
                ) * Fy

            self.f[i] += force_term

    def stream(self):
        """Simpler and clearer streaming with bounce-back"""
        # Temporary array for new distributions
        f_new = np.zeros_like(self.f)

        for i in range(Q):
            cx, cy = c[i, 0], c[i, 1]
            opp = opposite[i]

            for x in range(self.nx):
                for y in range(self.ny):
                    # Skip obstacle nodes (nothing streams from them)
                    if self.obstacle[x, y]:
                        continue

                    # Calculate destination
                    x_new = (x + cx) % self.nx
                    y_new = y + cy

                    # If destination is valid
                    if 0 <= y_new < self.ny:
                        if self.obstacle[x_new, y_new]:
                            # BOUNCE BACK at cylinder
                            f_new[opp, x, y] = self.f[i, x, y]
                        else:
                            # Normal streaming
                            f_new[i, x_new, y_new] = self.f[i, x, y]
                    else:
                        # Destination is at wall - wall BC will handle
                        # For now, treat as bounce-back (simple approach)
                        f_new[opp, x, y] = self.f[i, x, y]

        # Update distributions
        self.f = f_new

    def calculate_macroscopic(self):
        """Calculate macroscopic variables - OPTIMIZED"""
        # Calculate density (sum over all directions)
        self.rho = np.sum(self.f, axis=0)

        # IMPORTANT: Handle obstacle nodes separately
        fluid_mask = ~self.obstacle

        # Calculate momentum
        mom_x = np.sum(self.f * c[:, 0].reshape(Q, 1, 1), axis=0)
        mom_y = np.sum(self.f * c[:, 1].reshape(Q, 1, 1), axis=0)

        # Initialize velocities to zero
        self.u[0].fill(0.0)
        self.u[1].fill(0.0)

        # Calculate velocities only for fluid nodes
        if np.any(fluid_mask):
            rho_fluid = self.rho[fluid_mask]

            # Ensure no division by zero (use a small epsilon)
            rho_fluid = np.maximum(rho_fluid, 1e-10)

            # Calculate velocities with Guo correction
            self.u[0][fluid_mask] = (mom_x[fluid_mask] + 0.5 * body_force[0]) / rho_fluid
            self.u[1][fluid_mask] = (mom_y[fluid_mask] + 0.5 * body_force[1]) / rho_fluid

        # Explicitly set obstacle nodes
        self.u[0][self.obstacle] = 0.0
        self.u[1][self.obstacle] = 0.0
        self.rho[self.obstacle] = rho0  # Reset density at obstacles

    def collide(self):
        # Apply body force
        self.apply_body_force()

        # Recalculate macroscopic after force application
        self.calculate_macroscopic()

        # Calculate equilibrium distribution
        self.calculate_equilibrium()

        # Collision: relax toward equilibrium
        self.f -= omega * (self.f - self.f_eq)

        # Ensure distributions at obstacle nodes are zero
        for i in range(Q):
            self.f[i][self.obstacle] = 0.0

    def step(self):
        """LBM step"""
        self.stream()
        apply_boundary_conditions(self)
        self.collide()

    def get_velocity_magnitude(self):
        """Get velocity magnitude for visualization"""
        velocity_magnitude = np.sqrt(self.u[0] ** 2 + self.u[1] ** 2)
        velocity_magnitude[self.obstacle] = 0.0
        return velocity_magnitude
