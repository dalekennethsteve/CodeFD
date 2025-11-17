import matplotlib.pyplot as plt
import numpy as np


def plot_velocity_field(lb, step, save_path=None):
    """Plot velocity field"""
    plt.figure(figsize=(12, 4))

    # Create meshgrid
    X, Y = np.meshgrid(np.arange(lb.nx), np.arange(lb.ny), indexing='ij')

    # Plot velocity magnitude
    velocity_magnitude = np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2)

    plt.subplot(1, 2, 1)
    plt.contourf(X, Y, velocity_magnitude, levels=50, cmap='jet')
    plt.colorbar(label='Velocity Magnitude')
    plt.title(f'Velocity Magnitude - Step {step}')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')

    # Plot streamlines
    plt.subplot(1, 2, 2)
    plt.streamplot(X.T, Y.T, lb.u[0].T, lb.u[1].T, density=1, color='b', linewidth=0.5)
    plt.title(f'Streamlines - Step {step}')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')

    plt.tight_layout()

    if save_path:
        plt.savefig(f'{save_path}/step_{step:06d}.png', dpi=150, bbox_inches='tight')

    plt.show(block=False)
    plt.pause(0.1)


