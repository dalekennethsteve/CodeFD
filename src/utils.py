import numpy as np

from src.parameters import tau


def check_convergence(current_u, previous_u, tolerance=1e-6):

    # Calculate velocity magnitude for both time steps
    current_vel_mag = np.sqrt(current_u[0] ** 2 + current_u[1] ** 2)
    previous_vel_mag = np.sqrt(previous_u[0] ** 2 + previous_u[1] ** 2)

    # Calculate maximum change
    max_change = np.max(np.abs(current_vel_mag - previous_vel_mag))

    # Check if converged
    converged = max_change < tolerance

    return converged, max_change

