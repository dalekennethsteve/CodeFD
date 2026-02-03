import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle, Rectangle

from src.parameters import cylinder_x, cylinder_y, cylinder_radius

#from parameters import cylinder_x, cylinder_radius, cylinder_y


def plot_velocity_field(lb, step, save_path=None):
    """Plot velocity field"""
    plt.figure(figsize=(10, 5))

    # Create meshgrid
    X, Y = np.meshgrid(np.arange(lb.nx), np.arange(lb.ny), indexing='ij')

    # Plot velocity magnitude
    velocity_magnitude = np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2)

    plt.contourf(X, Y, velocity_magnitude, levels=50, cmap='jet')

    # Add cylinder obstacle
    cylinder = Circle((cylinder_x, cylinder_y), cylinder_radius,
                      linewidth=2, edgecolor='black', facecolor='none')
    plt.gca().add_patch(cylinder)

    # Commented out square obstacle for later reuse
    
    # Add square obstacle
    square_half_size = cylinder_radius
    square_x = cylinder_x - square_half_size
    square_y = cylinder_y - square_half_size
    square_width = 2 * square_half_size
    square = Rectangle((square_x, square_y), square_width, square_width,
                       linewidth=2, edgecolor='black', facecolor='none')
    plt.gca().add_patch(square)
    

    plt.colorbar(label='Velocity Magnitude')
    plt.title(f'Velocity Magnitude - Step {step}')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')

    plt.tight_layout()

    if save_path:
        plt.savefig(f'{save_path}/step_{step:06d}.png', dpi=150, bbox_inches='tight')

    plt.show(block=False)
    plt.pause(0.1)


def plot_pressure_field(lb, step, save_path=None):
    """Plot pressure field"""
    plt.figure(figsize=(10, 5))

    # Create meshgrid
    X, Y = np.meshgrid(np.arange(lb.nx), np.arange(lb.ny), indexing='ij')

    # Calculate pressure (from density)
    # In LBM, pressure = cs^2 * density, where cs^2 = 1/3
    pressure = (1 / 3) * lb.rho

    # Plot pressure field
    plt.contourf(X, Y, pressure, levels=50, cmap='jet')

    # Add cylinder obstacle
    cylinder = Circle((cylinder_x, cylinder_y), cylinder_radius,
                      linewidth=2, edgecolor='black', facecolor='none')
    plt.gca().add_patch(cylinder)

    # Commented out square obstacle for later reuse
    
    # Add square obstacle (same as velocity plot)
    square_half_size = cylinder_radius
    square_x = cylinder_x - square_half_size
    square_y = cylinder_y - square_half_size
    square_width = 2 * square_half_size
    square = Rectangle((square_x, square_y), square_width, square_width,
                       linewidth=2, edgecolor='black', facecolor='none')
    plt.gca().add_patch(square)
    

    plt.colorbar(label='Pressure')
    plt.title(f'Pressure Field - Step {step}')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')

    plt.tight_layout()

    if save_path:
        plt.savefig(f'{save_path}/pressure_step_{step:06d}.png', dpi=150, bbox_inches='tight')

    plt.show(block=False)
    plt.pause(0.1)


def plot_velocity_with_streamlines(lb, step, save_path=None):
    plt.figure(figsize=(10, 5))

    # Create meshgrid
    X, Y = np.meshgrid(np.arange(lb.nx), np.arange(lb.ny), indexing='ij')

    # Plot velocity magnitude as background
    velocity_magnitude = np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2)
    plt.contourf(X, Y, velocity_magnitude, levels=50, cmap='jet', alpha=0.7)

    # ===== SQUARE mask =====
    square_side = cylinder_radius * 2
    square_left = cylinder_x - cylinder_radius
    square_right = cylinder_x + cylinder_radius
    square_bottom = cylinder_y - cylinder_radius
    square_top = cylinder_y + cylinder_radius

    square_mask = ((X >= square_left) & (X <= square_right) &
                   (Y >= square_bottom) & (Y <= square_top))

    obstacle_mask = square_mask
    # ================================================

    # Create copies of velocities
    u_masked = lb.u[0].copy()
    v_masked = lb.u[1].copy()

    # Set velocities to NaN INSIDE SQUARE ONLY
    u_masked[obstacle_mask] = np.nan
    v_masked[obstacle_mask] = np.nan

    # ===== NO BUFFER ZONE =====

    # Downsample for streamlines
    skip = 3
    X_ds = X[::skip, ::skip]
    Y_ds = Y[::skip, ::skip]
    u_ds = u_masked[::skip, ::skip]
    v_ds = v_masked[::skip, ::skip]

    # Overlay streamlines - use BLACK for better contrast
    plt.streamplot(X_ds.T, Y_ds.T, u_ds.T, v_ds.T,
                   color='black',  # Black shows better than white on colored background
                   linewidth=0.8,  # Slightly thicker
                   density=2.0,  # Increased density
                   arrowsize=0.7,
                   broken_streamlines=False,
                   minlength=0.5)

    # Add SQUARE obstacle (semi-transparent)
    square = Rectangle((square_left, square_bottom),
                       square_side, square_side,
                       linewidth=2, edgecolor='black', facecolor='white', alpha=0.6)  # Semi-transparent white
    plt.gca().add_patch(square)

    plt.colorbar(label='Velocity Magnitude')
    plt.title(f'Velocity with Streamlines - Step {step}')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')

    plt.tight_layout()

    if save_path:
        plt.savefig(f'{save_path}/velocity_streamlines_step_{step:06d}.png', dpi=150, bbox_inches='tight')

    plt.show(block=False)
    plt.pause(0.1)



def plot_streamlines_zoom(lb, step, save_path=None):
    """
    Zoomed streamlines plot around a static square obstacle.
    Square side length = 2 * cylinder_radius.
    """
    plt.figure(figsize=(10, 5))

    # Coordinate arrays
    x = np.arange(lb.nx)
    y = np.arange(lb.ny)

    # Square mask
    X, Y = np.meshgrid(x, y, indexing='ij')
    R = cylinder_radius
    cx, cy = cylinder_x, cylinder_y

    square_mask = (np.abs(X - cx) <= R) & (np.abs(Y - cy) <= R)

    # Mask velocities inside square
    u_masked = lb.u[0].copy()
    v_masked = lb.u[1].copy()
    u_masked[square_mask] = np.nan
    v_masked[square_mask] = np.nan

    # Downsampling (keep dense for zoom)
    skip = 1
    x_ds = x[::skip]
    y_ds = y[::skip]
    u_ds = u_masked[::skip, ::skip]
    v_ds = v_masked[::skip, ::skip]

    # Streamlines
    plt.streamplot(
        x_ds, y_ds,
        u_ds.T, v_ds.T,
        color='darkblue',
        linewidth=0.6,
        density=3.0,
        arrowsize=0.6,
        broken_streamlines=False
    )

    # Draw square obstacle
    square = Rectangle(
        (cx - R, cy - R),
        2 * R, 2 * R,
        linewidth=2,
        edgecolor='red',
        facecolor='red',
        alpha=0.8
    )
    plt.gca().add_patch(square)

    # Zoom around square
    zoom_radius = 3 * R
    plt.xlim(cx - zoom_radius, cx + zoom_radius)
    plt.ylim(cy - zoom_radius, cy + zoom_radius)

    # Labels and title
    plt.title(f'Streamlines around a Square Obstacle – Step {step}')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save if requested
    if save_path:
        import os
        os.makedirs(save_path, exist_ok=True)
        filename = f'{save_path}/static_square_streamlines_zoom_step_{step:06d}.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"  Saved: {filename}")

    plt.show(block=False)
    plt.pause(0.1)
