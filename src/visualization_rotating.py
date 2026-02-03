"""
main_rotating_minimal.py - Minimal rotating cylinder with only velocity plots
"""
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from rotating_cylinder import RotatingCylinder
from parameters import *


def plot_velocity_rotating(lb, step, omega_rot=0.0, save_path=None):
    """
    Plot velocity field for rotating cylinder
    Fixed version - save_path should be a string folder path, not omega_rot
    """
    plt.figure(figsize=(10, 5))

    # Create meshgrid
    X, Y = np.meshgrid(np.arange(lb.nx), np.arange(lb.ny), indexing='ij')

    # Plot velocity magnitude
    velocity_magnitude = np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2)

    plt.contourf(X, Y, velocity_magnitude, levels=50, cmap='jet')

    # Add cylinder outline
    from matplotlib.patches import Circle
    cylinder = Circle((cylinder_x, cylinder_y), cylinder_radius,
                      linewidth=2, edgecolor='black', facecolor='none')
    plt.gca().add_patch(cylinder)

    plt.colorbar(label='Velocity Magnitude')

    # Title with rotation info
    if omega_rot == 0:
        plt.title(f'Velocity Magnitude - Static Cylinder - Step {step}')
    else:
        plt.title(f'Velocity Magnitude - Rotating Cylinder (ω={omega_rot:.3f}) - Step {step}')

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')

    plt.tight_layout()

    # FIXED: Check if save_path is provided and is a valid directory
    if save_path:
        import os
        # Create directory if it doesn't exist
        os.makedirs(save_path, exist_ok=True)

        # Generate filename
        if omega_rot == 0:
            filename = f'{save_path}/static_step_{step:06d}.png'
        else:
            # Use absolute value for filename to avoid negative signs in folder names
            omega_str = f'{abs(omega_rot):.3f}'.replace('.', 'p')
            direction = 'ccw' if omega_rot > 0 else 'cw'
            filename = f'{save_path}/rotating_{omega_str}_{direction}_step_{step:06d}.png'

        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"  Saved: {filename}")

    plt.show(block=False)
    plt.pause(0.1)


def plot_pressure_rotating(lb, step, omega_rot=0.0, save_path=None):
    """
    Plot pressure field for rotating cylinder
    Same layout as velocity plot
    """
    plt.figure(figsize=(10, 5))

    # Create meshgrid
    X, Y = np.meshgrid(np.arange(lb.nx), np.arange(lb.ny), indexing='ij')

    # Calculate pressure (from density)
    # In LBM, pressure = cs^2 * density, where cs^2 = 1/3
    pressure = (1 / 3) * lb.rho

    # Plot pressure field with diverging colormap (better for pressure)
    plt.contourf(X, Y, pressure, levels=50, cmap='jet')

    # Add cylinder outline
    from matplotlib.patches import Circle
    cylinder = Circle((cylinder_x, cylinder_y), cylinder_radius,
                      linewidth=2, edgecolor='black', facecolor='none')
    plt.gca().add_patch(cylinder)

    plt.colorbar(label='Pressure')

    # Title with rotation info
    if omega_rot == 0:
        plt.title(f'Pressure Field - Static Cylinder - Step {step}')
    else:
        plt.title(f'Pressure Field - Rotating Cylinder (ω={omega_rot:.3f}) - Step {step}')

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')

    plt.tight_layout()

    # Save if save_path is provided
    if save_path:
        import os
        # Create directory if it doesn't exist
        os.makedirs(save_path, exist_ok=True)

        # Generate filename
        if omega_rot == 0:
            filename = f'{save_path}/static_pressure_step_{step:06d}.png'
        else:
            # Use absolute value for filename to avoid negative signs in folder names
            omega_str = f'{abs(omega_rot):.3f}'.replace('.', 'p')
            direction = 'ccw' if omega_rot > 0 else 'cw'
            filename = f'{save_path}/rotating_pressure_{omega_str}_{direction}_step_{step:06d}.png'

        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"  Saved: {filename}")

    plt.show(block=False)
    plt.pause(0.1)



def plot_combined_rotating(lb, step, omega_rot=0.0, save_path=None):
    """
    Combined plot showing velocity magnitude and streamlines
    """
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))

    # Create meshgrid
    X, Y = np.meshgrid(np.arange(lb.nx), np.arange(lb.ny), indexing='ij')

    # Left subplot: Velocity magnitude
    velocity_magnitude = np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2)
    im1 = axes[0].contourf(X, Y, velocity_magnitude, levels=50, cmap='jet')

    # Add cylinder to velocity plot
    cylinder1 = Circle((cylinder_x, cylinder_y), cylinder_radius,
                       linewidth=2, edgecolor='black', facecolor='none')
    axes[0].add_patch(cylinder1)

    plt.colorbar(im1, ax=axes[0], label='Velocity Magnitude')
    axes[0].set_title('Velocity Magnitude')
    axes[0].set_xlabel('X')
    axes[0].set_ylabel('Y')
    axes[0].axis('equal')

    # Right subplot: Streamlines
    axes[1].streamplot(X.T, Y.T, lb.u[0].T, lb.u[1].T,
                       color='darkblue', linewidth=0.8, density=1.5, arrowsize=0.8)

    # Add cylinder to streamline plot
    cylinder2 = Circle((cylinder_x, cylinder_y), cylinder_radius,
                       linewidth=2, edgecolor='black', facecolor='black', alpha=0.8)
    axes[1].add_patch(cylinder2)

    axes[1].set_title('Streamlines')
    axes[1].set_xlabel('X')
    axes[1].set_ylabel('Y')
    axes[1].axis('equal')

    # Main title
    if omega_rot == 0:
        fig.suptitle(f'Flow Field - Static Cylinder - Step {step}', fontsize=14)
    else:
        fig.suptitle(f'Flow Field - Rotating Cylinder (ω={omega_rot:.3f}) - Step {step}', fontsize=14)

    plt.tight_layout()

    # Save if save_path is provided
    if save_path:
        import os
        os.makedirs(save_path, exist_ok=True)

        if omega_rot == 0:
            filename = f'{save_path}/static_combined_step_{step:06d}.png'
        else:
            omega_str = f'{abs(omega_rot):.3f}'.replace('.', 'p')
            direction = 'ccw' if omega_rot > 0 else 'cw'
            filename = f'{save_path}/rotating_combined_{omega_str}_{direction}_step_{step:06d}.png'

        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"  Saved combined plot: {filename}")

    plt.show(block=False)
    plt.pause(0.1)


def plot_force_history(force_history, omega_rot, save_path=None):
    """
    Plot drag and lift forces over time
    force_history: list of (step, drag, lift) tuples
    """
    if not force_history:
        print("No force history to plot")
        return

    steps = [fh[0] for fh in force_history]
    drags = [fh[1] for fh in force_history]
    lifts = [fh[2] for fh in force_history]

    plt.figure(figsize=(12, 5))

    # Plot drag force
    plt.subplot(1, 2, 1)
    plt.plot(steps, drags, 'r-', linewidth=2, label='Drag Force')
    plt.xlabel('Time Step')
    plt.ylabel('Drag Force (lattice units)')
    plt.title(f'Drag Force Evolution (ω={omega_rot:.3f})')
    plt.grid(True, alpha=0.3)
    plt.legend()

    # Plot lift force
    plt.subplot(1, 2, 2)
    plt.plot(steps, lifts, 'b-', linewidth=2, label='Lift Force')
    plt.xlabel('Time Step')
    plt.ylabel('Lift Force (lattice units)')
    plt.title(f'Lift Force Evolution (ω={omega_rot:.3f})')
    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.tight_layout()

    if save_path:
        import os
        os.makedirs(save_path, exist_ok=True)
        filename = f'{save_path}/forces_omega_{abs(omega_rot):.3f}.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"  Saved force plot: {filename}")

    plt.show(block=False)
    plt.pause(0.1)



def plot_combined_rotating(lb, step, omega_rot=0.0, save_path=None):
    """
    Combined plot showing velocity magnitude and streamlines
    """
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))

    # Create meshgrid
    X, Y = np.meshgrid(np.arange(lb.nx), np.arange(lb.ny), indexing='ij')

    # Left subplot: Velocity magnitude
    velocity_magnitude = np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2)
    im1 = axes[0].contourf(X, Y, velocity_magnitude, levels=50, cmap='jet')

    # Add cylinder to velocity plot
    cylinder1 = Circle((cylinder_x, cylinder_y), cylinder_radius,
                       linewidth=2, edgecolor='black', facecolor='none')
    axes[0].add_patch(cylinder1)

    plt.colorbar(im1, ax=axes[0], label='Velocity Magnitude')
    axes[0].set_title('Velocity Magnitude')
    axes[0].set_xlabel('X')
    axes[0].set_ylabel('Y')
    axes[0].axis('equal')

    # Right subplot: Streamlines
    axes[1].streamplot(X.T, Y.T, lb.u[0].T, lb.u[1].T,
                       color='darkblue', linewidth=0.8, density=1.5, arrowsize=0.8)

    # Add cylinder to streamline plot
    cylinder2 = Circle((cylinder_x, cylinder_y), cylinder_radius,
                       linewidth=2, edgecolor='black', facecolor='black', alpha=0.8)
    axes[1].add_patch(cylinder2)

    axes[1].set_title('Streamlines')
    axes[1].set_xlabel('X')
    axes[1].set_ylabel('Y')
    axes[1].axis('equal')

    # Main title - FIXED: positive omega = CW
    if omega_rot == 0:
        fig.suptitle(f'Flow Field - Static Cylinder - Step {step}', fontsize=14)
    else:
        rotation_dir = "CW" if omega_rot > 0 else "CCW"  # FIXED
        fig.suptitle(f'Flow Field - Rotating Cylinder (ω={omega_rot:.3f}, {rotation_dir}) - Step {step}', fontsize=14)

    plt.tight_layout()

    # Save if save_path is provided
    if save_path:
        import os
        os.makedirs(save_path, exist_ok=True)

        if omega_rot == 0:
            filename = f'{save_path}/static_combined_step_{step:06d}.png'
        else:
            omega_str = f'{abs(omega_rot):.3f}'.replace('.', 'p')
            direction = 'cw' if omega_rot > 0 else 'ccw'  # FIXED
            filename = f'{save_path}/rotating_combined_{omega_str}_{direction}_step_{step:06d}.png'

        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"  Saved combined plot: {filename}")

    plt.show(block=False)
    plt.pause(0.1)


# FIXED streamlines functions below:

def plot_streamlines_rotating(lb, step, omega_rot=0.0, save_path=None):
    """
    Plot streamlines for rotating cylinder
    Same layout as velocity plot
    """
    plt.figure(figsize=(10, 5))

    # Create coordinate arrays
    x = np.arange(lb.nx)
    y = np.arange(lb.ny)

    # Create cylinder mask to hide streamlines inside
    X, Y = np.meshgrid(x, y, indexing='ij')
    dx = X - cylinder_x
    dy = Y - cylinder_y
    cylinder_mask = (dx ** 2 + dy ** 2) <= cylinder_radius ** 2

    # Mask velocities inside cylinder
    u_masked = lb.u[0].copy()
    v_masked = lb.u[1].copy()
    u_masked[cylinder_mask] = np.nan
    v_masked[cylinder_mask] = np.nan

    # Add small buffer
    buffer_radius = cylinder_radius + 1
    buffer_mask = (dx ** 2 + dy ** 2) <= buffer_radius ** 2
    u_masked[buffer_mask] = np.nan
    v_masked[buffer_mask] = np.nan

    # Downsample for cleaner streamlines
    skip = 2
    x_ds = x[::skip]
    y_ds = y[::skip]
    u_ds = u_masked[::skip, ::skip]
    v_ds = v_masked[::skip, ::skip]

    # Plot streamlines
    plt.streamplot(x_ds, y_ds,
                   u_ds.T,
                   v_ds.T,
                   color='darkblue',
                   linewidth=0.8,
                   density=2.0,
                   arrowsize=0.8,
                   broken_streamlines=False)

    # Add cylinder outline
    cylinder = Circle((cylinder_x, cylinder_y), cylinder_radius,
                      linewidth=2, edgecolor='black', facecolor='black', alpha=0.9)
    plt.gca().add_patch(cylinder)

    # Title with rotation info - FIXED: positive omega = CW
    if omega_rot == 0:
        plt.title(f'Streamlines - Static Cylinder - Step {step}')
    else:
        rotation_dir = "CW" if omega_rot > 0 else "CCW"  # FIXED
        plt.title(f'Streamlines - Rotating Cylinder (ω={omega_rot:.3f}, {rotation_dir}) - Step {step}')

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')

    plt.tight_layout()

    # Save if save_path is provided - FIXED: positive omega = cw
    if save_path:
        import os
        os.makedirs(save_path, exist_ok=True)

        if omega_rot == 0:
            filename = f'{save_path}/static_streamlines_step_{step:06d}.png'
        else:
            omega_str = f'{abs(omega_rot):.3f}'.replace('.', 'p')
            direction = 'cw' if omega_rot > 0 else 'ccw'  # FIXED
            filename = f'{save_path}/rotating_streamlines_{omega_str}_{direction}_step_{step:06d}.png'

        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"  Saved: {filename}")

    plt.show(block=False)
    plt.pause(0.1)


def plot_streamlines_with_velocity_color(lb, step, omega_rot=0.0, save_path=None):
    """
    Plot streamlines colored by velocity magnitude for rotating cylinder
    """
    plt.figure(figsize=(10, 5))

    # Create coordinate arrays
    x = np.arange(lb.nx)
    y = np.arange(lb.ny)

    # Calculate velocity magnitude
    velocity_magnitude = np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2)

    # Create cylinder mask
    X, Y = np.meshgrid(x, y, indexing='ij')
    dx = X - cylinder_x
    dy = Y - cylinder_y
    cylinder_mask = (dx ** 2 + dy ** 2) <= cylinder_radius ** 2

    # Mask velocities inside cylinder
    u_masked = lb.u[0].copy()
    v_masked = lb.u[1].copy()
    vel_masked = velocity_magnitude.copy()

    u_masked[cylinder_mask] = np.nan
    v_masked[cylinder_mask] = np.nan
    vel_masked[cylinder_mask] = np.nan

    # Downsample
    skip = 2
    x_ds = x[::skip]
    y_ds = y[::skip]
    u_ds = u_masked[::skip, ::skip]
    v_ds = v_masked[::skip, ::skip]
    vel_ds = vel_masked[::skip, ::skip]

    # Plot streamlines colored by velocity
    strm = plt.streamplot(x_ds, y_ds,
                          u_ds.T,
                          v_ds.T,
                          color=vel_ds.T,
                          cmap='jet',
                          linewidth=1.0,
                          density=2.0,
                          arrowsize=0.8,
                          broken_streamlines=False)

    # Add cylinder
    cylinder = Circle((cylinder_x, cylinder_y), cylinder_radius,
                      linewidth=2, edgecolor='black', facecolor='black', alpha=0.9)
    plt.gca().add_patch(cylinder)

    # Add colorbar
    cbar = plt.colorbar(strm.lines, label='Velocity Magnitude')

    # Title - FIXED: positive omega = CW
    if omega_rot == 0:
        plt.title(f'Streamlines with Velocity Color - Static Cylinder - Step {step}')
    else:
        rotation_dir = "CW" if omega_rot > 0 else "CCW"  # FIXED
        plt.title(
            f'Streamlines with Velocity Color - Rotating Cylinder (ω={omega_rot:.3f}, {rotation_dir}) - Step {step}')

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')

    plt.tight_layout()

    # Save if save_path is provided - FIXED: positive omega = cw
    if save_path:
        import os
        os.makedirs(save_path, exist_ok=True)

        if omega_rot == 0:
            filename = f'{save_path}/static_streamlines_color_step_{step:06d}.png'
        else:
            omega_str = f'{abs(omega_rot):.3f}'.replace('.', 'p')
            direction = 'cw' if omega_rot > 0 else 'ccw'  # FIXED
            filename = f'{save_path}/rotating_streamlines_color_{omega_str}_{direction}_step_{step:06d}.png'

        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"  Saved: {filename}")

    plt.show(block=False)
    plt.pause(0.1)


def plot_streamlines_zoom_rotating(lb, step, omega_rot=0.0, save_path=None):
    """
    Zoomed streamlines plot around rotating cylinder
    """
    plt.figure(figsize=(10, 5))

    # Create coordinate arrays
    x = np.arange(lb.nx)
    y = np.arange(lb.ny)

    # Create cylinder mask
    X, Y = np.meshgrid(x, y, indexing='ij')
    dx = X - cylinder_x
    dy = Y - cylinder_y
    cylinder_mask = (dx ** 2 + dy ** 2) <= cylinder_radius ** 2

    # Mask velocities
    u_masked = lb.u[0].copy()
    v_masked = lb.u[1].copy()
    u_masked[cylinder_mask] = np.nan
    v_masked[cylinder_mask] = np.nan

    # Downsample less for zoomed view
    skip = 1
    x_ds = x[::skip]
    y_ds = y[::skip]
    u_ds = u_masked[::skip, ::skip]
    v_ds = v_masked[::skip, ::skip]

    # Plot streamlines
    plt.streamplot(x_ds, y_ds,
                   u_ds.T,
                   v_ds.T,
                   color='darkblue',
                   linewidth=0.6,
                   density=3.0,  # Higher density for zoomed view
                   arrowsize=0.6,
                   broken_streamlines=False)

    # Add cylinder
    cylinder = Circle((cylinder_x, cylinder_y), cylinder_radius,
                      linewidth=2, edgecolor='red', facecolor='red', alpha=0.8)
    plt.gca().add_patch(cylinder)

    # Zoom around cylinder
    zoom_radius = cylinder_radius * 3
    plt.xlim(cylinder_x - zoom_radius, cylinder_x + zoom_radius)
    plt.ylim(cylinder_y - zoom_radius, cylinder_y + zoom_radius)

    # Title - FIXED: positive omega = CW
    if omega_rot == 0:
        plt.title(f'Zoomed Streamlines - Static Cylinder - Step {step}')
    else:
        rotation_dir = "CW" if omega_rot > 0 else "CCW"  # FIXED
        plt.title(f'Zoomed Streamlines - Rotating Cylinder (ω={omega_rot:.3f}, {rotation_dir}) - Step {step}')

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    plt.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save if save_path is provided - FIXED: positive omega = cw
    if save_path:
        import os
        os.makedirs(save_path, exist_ok=True)

        if omega_rot == 0:
            filename = f'{save_path}/static_streamlines_zoom_step_{step:06d}.png'
        else:
            omega_str = f'{abs(omega_rot):.3f}'.replace('.', 'p')
            direction = 'cw' if omega_rot > 0 else 'ccw'  # FIXED
            filename = f'{save_path}/rotating_streamlines_zoom_{omega_str}_{direction}_step_{step:06d}.png'

        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"  Saved: {filename}")

    plt.show(block=False)
    plt.pause(0.1)
