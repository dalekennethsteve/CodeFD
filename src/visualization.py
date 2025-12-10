# visualization.py - Updated to show cylinder
import matplotlib.pyplot as plt
import numpy as np

from src.parameters import cylinder_x, cylinder_radius, cylinder_y


def plot_velocity_field(lb, step, save_path=None):
    """Plot velocity field"""
    plt.figure(figsize=(10, 5))

    # Create meshgrid
    X, Y = np.meshgrid(np.arange(lb.nx), np.arange(lb.ny), indexing='ij')

    # Plot velocity magnitude
    velocity_magnitude = np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2)

    plt.contourf(X, Y, velocity_magnitude, levels=50, cmap='jet')

    # Add cylinder outline (ADDED)
    theta = np.linspace(0, 2 * np.pi, 100)
    circle_x = cylinder_x + cylinder_radius * np.cos(theta)
    circle_y = cylinder_y + cylinder_radius * np.sin(theta)
    plt.plot(circle_x, circle_y, 'k-', linewidth=2, label='Cylinder')

    plt.colorbar(label='Velocity Magnitude')
    plt.title(f'Velocity Magnitude - Step {step}')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    plt.legend()

    plt.tight_layout()

    if save_path:
        plt.savefig(f'{save_path}/step_{step:06d}.png', dpi=150, bbox_inches='tight')

    plt.show(block=False)
    plt.pause(0.1)