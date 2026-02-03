"""
debug_rotation.py - Debug why rotation isn't working
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from rotating_cylinder import RotatingCylinder
from lattice_boltzmann import LatticeBoltzmann
from parameters import *


def check_rotation_implementation():
    print("DEBUGGING ROTATION IMPLEMENTATION")
    print("=" * 60)

    # Test with stronger rotation
    omega_rot = 0.1

    print(f"Testing ω = {omega_rot}")

    # Create rotating cylinder
    lb = RotatingCylinder(omega_rot=omega_rot)

    print(f"\n1. Checking boundary_velocity values:")
    print(f"   Shape: {lb.boundary_velocity.shape}")
    print(f"   Max velocity: {np.max(np.abs(lb.boundary_velocity)):.6f}")
    print(f"   Min velocity: {np.min(np.abs(lb.boundary_velocity)):.6f}")

    # Check velocities at specific points on cylinder surface
    print(f"\n2. Velocities at cylinder surface points:")

    # Top of cylinder
    x, y = cylinder_x, cylinder_y + cylinder_radius
    print(f"   Top ({x}, {y}): u={lb.boundary_velocity[0, x, y]:.4f}, v={lb.boundary_velocity[1, x, y]:.4f}")

    # Right of cylinder
    x, y = cylinder_x + cylinder_radius, cylinder_y
    print(f"   Right ({x}, {y}): u={lb.boundary_velocity[0, x, y]:.4f}, v={lb.boundary_velocity[1, x, y]:.4f}")

    # Bottom of cylinder
    x, y = cylinder_x, cylinder_y - cylinder_radius
    print(f"   Bottom ({x}, {y}): u={lb.boundary_velocity[0, x, y]:.4f}, v={lb.boundary_velocity[1, x, y]:.4f}")

    # Left of cylinder
    x, y = cylinder_x - cylinder_radius, cylinder_y
    print(f"   Left ({x}, {y}): u={lb.boundary_velocity[0, x, y]:.4f}, v={lb.boundary_velocity[1, x, y]:.4f}")

    # Run a few steps and check if velocities change
    print(f"\n3. Running 10 steps and checking fluid velocities near cylinder:")

    for step in range(10):
        lb.step()

        if step % 2 == 0:
            # Check fluid velocities just outside cylinder
            x_right = cylinder_x + cylinder_radius + 1
            x_left = cylinder_x - cylinder_radius - 1
            y_center = cylinder_y

            u_right = lb.u[0, x_right, y_center]
            u_left = lb.u[0, x_left, y_center]

            print(f"   Step {step}: u_right={u_right:.6f}, u_left={u_left:.6f}, diff={u_right - u_left:.6f}")

    print(f"\n4. Visualizing boundary velocities:")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Plot boundary velocity magnitude
    ax1 = axes[0]
    boundary_vel_mag = np.sqrt(lb.boundary_velocity[0] ** 2 + lb.boundary_velocity[1] ** 2)
    im1 = ax1.imshow(boundary_vel_mag.T, origin='lower', cmap='hot', aspect='auto')

    # Add cylinder outline
    cylinder1 = Circle((cylinder_x, cylinder_y), cylinder_radius,
                       color='cyan', fill=False, linewidth=2)
    ax1.add_patch(cylinder1)
    ax1.set_title('Boundary Velocity Magnitude on Cylinder')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    plt.colorbar(im1, ax=ax1, label='Velocity')

    # Plot boundary velocity direction (quiver)
    ax2 = axes[1]

    # Downsample for clarity
    skip = 2
    X, Y = np.meshgrid(np.arange(0, nx, skip), np.arange(0, ny, skip), indexing='ij')
    U = lb.boundary_velocity[0, ::skip, ::skip]
    V = lb.boundary_velocity[1, ::skip, ::skip]

    # Create a mask for cylinder region
    mask = np.zeros_like(U, dtype=bool)
    for i in range(U.shape[0]):
        for j in range(U.shape[1]):
            x = i * skip
            y = j * skip
            dx = x - cylinder_x
            dy = y - cylinder_y
            if dx * dx + dy * dy <= cylinder_radius ** 2:
                mask[i, j] = True

    # Only plot inside cylinder
    U[~mask] = 0
    V[~mask] = 0

    ax2.quiver(X, Y, U, V, scale=5, color='red')

    cylinder2 = Circle((cylinder_x, cylinder_y), cylinder_radius,
                       color='blue', fill=False, linewidth=2)
    ax2.add_patch(cylinder2)
    ax2.set_title('Boundary Velocity Direction (should be tangential)')
    ax2.set_xlabel('X')
    ax2.set_aspect('equal')

    plt.tight_layout()
    plt.show()

    print(f"\n5. What to look for:")
    print(f"   - Boundary velocities should be TANGENTIAL to cylinder")
    print(f"   - Top should have velocity to the RIGHT (for CCW rotation)")
    print(f"   - Right should have velocity DOWN")
    print(f"   - Bottom should have velocity LEFT")
    print(f"   - Left should have velocity UP")


if __name__ == "__main__":
    check_rotation_implementation()