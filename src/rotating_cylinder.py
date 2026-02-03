"""
rotating_cylinder_final.py - FINAL working version
"""
import numpy as np
from parameters import nx, ny, rho0, omega, body_force, cylinder_x, cylinder_y, cylinder_radius
from parameters import Q, c, w, opposite



class RotatingCylinder:
    def __init__(self, omega_rot=0.1):
        self.nx = nx
        self.ny = ny
        self.f = np.zeros((Q, nx, ny))
        self.f_eq = np.zeros((Q, nx, ny))
        self.rho = np.ones((nx, ny)) * rho0
        self.u = np.zeros((2, nx, ny))

        # Rotational parameter
        self.omega_rot = omega_rot

        # CRITICAL: Initialize boundary_velocity BEFORE create_cylinder_mask
        self.boundary_velocity = np.zeros((2, nx, ny))

        # Create obstacle mask
        self.obstacle = np.zeros((nx, ny), dtype=bool)
        self.create_cylinder_mask()

        self.initialize()

        max_vel = np.max(np.abs(self.boundary_velocity))
        print(f"✓ Rotating cylinder: ω = {omega_rot}")
        print(f"  Max boundary velocity: {max_vel:.4f}")
        print(f"  Boundary velocity shape: {self.boundary_velocity.shape}")

    def create_cylinder_mask(self):
        """Create mask with tangential velocities - MUST use self.boundary_velocity"""
        for x in range(self.nx):
            for y in range(self.ny):
                dx = x - cylinder_x
                dy = y - cylinder_y
                dist_sq = dx * dx + dy * dy

                if dist_sq <= cylinder_radius ** 2:
                    self.obstacle[x, y] = True

                    # Calculate tangential velocity
                    if dist_sq > 0:
                        distance = np.sqrt(dist_sq)
                        # Tangential velocity for CCW rotation
                        u_tangent = -self.omega_rot * dy / distance
                        v_tangent = self.omega_rot * dx / distance

                        # Store in boundary_velocity
                        self.boundary_velocity[0, x, y] = u_tangent
                        self.boundary_velocity[1, x, y] = v_tangent

    def initialize(self):
        """Initialize"""
        self.rho[:, :] = rho0
        self.u[0, :, :] = 0.0
        self.u[1, :, :] = 0.0
        self.calculate_equilibrium()
        self.f[:, :, :] = self.f_eq[:, :, :]

    def calculate_equilibrium(self):
        """Calculate equilibrium"""
        u_sqr = self.u[0] ** 2 + self.u[1] ** 2
        for i in range(Q):
            cu = c[i, 0] * self.u[0] + c[i, 1] * self.u[1]
            self.f_eq[i] = w[i] * self.rho * (1 + 3 * cu + 4.5 * cu ** 2 - 1.5 * u_sqr)

    def apply_body_force(self):
        """Apply body force"""
        Fx, Fy = body_force

        for i in range(Q):
            cu = c[i, 0] * self.u[0] + c[i, 1] * self.u[1]

            force_term = w[i] * (1 - 0.5 * omega) * (
                    3.0 * (c[i, 0] - self.u[0]) + 9.0 * cu * c[i, 0]
            ) * Fx

            if Fy != 0:
                force_term += w[i] * (1 - 0.5 * omega) * (
                        3.0 * (c[i, 1] - self.u[1]) + 9.0 * cu * c[i, 1]
                ) * Fy

            self.f[i] += force_term

    def stream(self):
        """Streaming with moving boundary"""
        f_new = np.zeros_like(self.f)

        for i in range(Q):
            cx, cy = c[i, 0], c[i, 1]
            opp = opposite[i]

            for x in range(self.nx):
                for y in range(self.ny):
                    if self.obstacle[x, y]:
                        continue

                    x_new = (x + cx) % self.nx
                    y_new = y + cy

                    if 0 <= y_new < self.ny:
                        if self.obstacle[x_new, y_new]:
                            # MOVING BOUNDARY - KEY FOR ROTATION
                            rho_local = self.rho[x, y]
                            u_wall = self.boundary_velocity[0, x_new, y_new]
                            v_wall = self.boundary_velocity[1, x_new, y_new]

                            # Simple moving boundary formula
                            cu_wall = c[opp, 0] * u_wall + c[opp, 1] * v_wall
                            f_new[opp, x, y] = self.f[i, x, y] - 2 * w[opp] * rho_local * (3 * cu_wall)
                        else:
                            f_new[i, x_new, y_new] = self.f[i, x, y]
                    else:
                        f_new[opp, x, y] = self.f[i, x, y]

        self.f = f_new

    def calculate_macroscopic(self):
        """Calculate macroscopic"""
        self.rho = np.sum(self.f, axis=0)

        fluid_mask = ~self.obstacle

        mom_x = np.sum(self.f * c[:, 0].reshape(Q, 1, 1), axis=0)
        mom_y = np.sum(self.f * c[:, 1].reshape(Q, 1, 1), axis=0)

        self.u[0].fill(0.0)
        self.u[1].fill(0.0)

        if np.any(fluid_mask):
            rho_fluid = self.rho[fluid_mask]
            rho_fluid = np.maximum(rho_fluid, 1e-10)

            self.u[0][fluid_mask] = (mom_x[fluid_mask] + 0.5 * body_force[0]) / rho_fluid
            self.u[1][fluid_mask] = (mom_y[fluid_mask] + 0.5 * body_force[1]) / rho_fluid

        # Set obstacle velocities
        self.u[0][self.obstacle] = self.boundary_velocity[0][self.obstacle]
        self.u[1][self.obstacle] = self.boundary_velocity[1][self.obstacle]
        self.rho[self.obstacle] = rho0

    def collide(self):
        """Collision"""
        self.apply_body_force()
        self.calculate_macroscopic()
        self.calculate_equilibrium()

        self.f -= omega * (self.f - self.f_eq)

        for i in range(Q):
            self.f[i][self.obstacle] = 0.0

    def step(self):
        """Complete step"""
        self.stream()
        self.apply_boundary_conditions()
        self.collide()

    def apply_boundary_conditions(self):
        """Boundary conditions"""
        # Periodic x
        for i in range(Q):
            self.f[i, -1, :] = self.f[i, 0, :]
            self.f[i, 0, :] = self.f[i, -1, :]

        # Bounce-back walls
        y = 0
        self.f[2, :, y] = self.f[4, :, y]
        self.f[5, :, y] = self.f[7, :, y]
        self.f[6, :, y] = self.f[8, :, y]

        y = self.ny - 1
        self.f[4, :, y] = self.f[2, :, y]
        self.f[7, :, y] = self.f[5, :, y]
        self.f[8, :, y] = self.f[6, :, y]

    def get_velocity_magnitude(self):
        """Get velocity"""
        velocity_magnitude = np.sqrt(self.u[0] ** 2 + self.u[1] ** 2)
        velocity_magnitude[self.obstacle] = 0.0
        return velocity_magnitude

    def calculate_forces(self):
        """
        Calculate lift and drag forces on the cylinder using momentum exchange method
        Returns: (drag_force, lift_force) in lattice units
        """
        drag_force = 0.0
        lift_force = 0.0

        # Scan all boundary nodes
        for x in range(self.nx):
            for y in range(self.ny):
                if not self.obstacle[x, y]:
                    continue

                # Check all neighboring fluid nodes
                for i in range(Q):
                    cx, cy = c[i, 0], c[i, 1]
                    x_fluid = (x - cx) % self.nx  # Fluid node streaming INTO obstacle
                    y_fluid = y - cy

                    # Check if neighbor is in bounds and is fluid
                    if 0 <= y_fluid < self.ny and not self.obstacle[x_fluid, y_fluid]:
                        # Calculate momentum exchange
                        opp = opposite[i]

                        # Post-collision distribution function at fluid node
                        f_i = self.f[i, x_fluid, y_fluid]
                        f_opp_eq = self.f_eq[opp, x_fluid, y_fluid]

                        # Bounce-back with moving wall contribution
                        delta_momentum = 2 * (f_i - 2 * w[opp] * self.rho[x_fluid, y_fluid] *
                                              (3 * (c[opp, 0] * self.boundary_velocity[0, x, y] +
                                                    c[opp, 1] * self.boundary_velocity[1, x, y])))

                        # Add to forces (force = momentum change per time step)
                        drag_force += delta_momentum * cx  # Force in x-direction
                        lift_force += delta_momentum * cy  # Force in y-direction

        return drag_force, lift_force

    def calculate_pressure_field(self):
        """
        Calculate pressure field from density
        Pressure = density * cs^2, where cs^2 = 1/3 for D2Q9
        """
        cs2 = 1.0 / 3.0  # Speed of sound squared in lattice units
        pressure = self.rho * cs2
        pressure[self.obstacle] = 0.0  # Set obstacle pressure to 0 for visualization
        return pressure

    def calculate_vorticity_field(self):
        """
        Calculate vorticity field: ω = ∂v/∂x - ∂u/∂y
        """
        vorticity = np.zeros((self.nx, self.ny))

        # Central difference (simple approximation)
        for x in range(1, self.nx - 1):
            for y in range(1, self.ny - 1):
                if not self.obstacle[x, y]:
                    dudy = (self.u[0, x, y + 1] - self.u[0, x, y - 1]) / 2.0
                    dvdx = (self.u[1, x + 1, y] - self.u[1, x - 1, y]) / 2.0
                    vorticity[x, y] = dvdx - dudy

        vorticity[self.obstacle] = 0.0
        return vorticity