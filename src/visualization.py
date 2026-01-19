# visualization.py - multi-cylinder version
import matplotlib.pyplot as plt
import numpy as np

from src.parameters import cylinders, cylinder_radius


def plot_velocity_field(lb, step, save_path=None):
    """Plot velocity field with multiple cylinders"""

    plt.figure(figsize=(10, 5))

    # Create meshgrid
    X, Y = np.meshgrid(
        np.arange(lb.nx),
        np.arange(lb.ny),
        indexing='ij'
    )

    # Velocity magnitude
    velocity_magnitude = np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2)

    plt.contourf(
        X, Y,
        velocity_magnitude,
        levels=50,
        cmap='jet'
    )

    # Plot ALL cylinder outlines
    theta = np.linspace(0, 2 * np.pi, 200)
    for i, (cx, cy) in enumerate(cylinders):
        circle_x = cx + cylinder_radius * np.cos(theta)
        circle_y = cy + cylinder_radius * np.sin(theta)

        plt.plot(
            circle_x,
            circle_y,
            'k-',
            linewidth=2,
            label='Cylinder' if i == 0 else None
        )

    plt.colorbar(label='Velocity Magnitude')
    plt.title(f'Velocity Magnitude - Step {step}')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')

    if cylinders:
        plt.legend()

    plt.tight_layout()

    if save_path:
        plt.savefig(
            f'{save_path}/step_{step:06d}.png',
            dpi=150,
            bbox_inches='tight'
        )

    plt.show(block=False)
    plt.pause(0.1)
