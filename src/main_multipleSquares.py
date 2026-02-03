# main.py
import numpy as np
import time

import lattice_boltzmann
from visualization import plot_velocity_field, plot_streamlines_zoom, plot_pressure_field
from utils import check_convergence
from parameters import cylinders, cylinder_type, cylinder_width, Re


def main():
    print("=" * 60)
    print("Initializing LBM simulation: Flow past multiple squares")
    print("=" * 60)

    print(f"Number of squares: {len(cylinders)}")
    print(f"Square width     : {cylinder_width}")

    for i, (cx, cy) in enumerate(cylinders, start=1):
        print(f"  Square {i}: center = ({cx}, {cy})")

    # Initialize solver
    lb = lattice_boltzmann()

    obstacle_nodes = np.sum(lb.obstacle)
    print(f"Total obstacle nodes: {obstacle_nodes}")

    # -------------------------------
    # Simulation control parameters
    # -------------------------------
    total_steps = 10000
    plot_interval = 2000
    convergence_check_interval = 200
    convergence_tolerance = 1e-6

    print(f"Running simulation for {total_steps} steps...")
    print(f"Reynolds number: {Re:.2f}")

    start_time = time.time()
    converged = False
    final_step = total_steps
    max_change = 0.0

    # -------------------------------
    # Main time loop
    # -------------------------------
    for step in range(total_steps):

        # Store previous velocity for convergence check
        if step % convergence_check_interval == 0:
            prev_u = lb.u.copy()

        lb.step()

        # Check convergence
        if step % convergence_check_interval == 0 and step > 0:
            converged, max_change = check_convergence(
                lb.u, prev_u, convergence_tolerance
            )
            if converged:
                final_step = step
                print(f"\nConvergence achieved at step {step}")
                print(f"Maximum velocity change: {max_change:.2e}")
                break

        # Runtime monitoring
        if step % plot_interval == 0:
            vel_mag = np.sqrt(lb.u[0]**2 + lb.u[1]**2)
            print(f"Step {step}/{total_steps}")
            print(f"  Max velocity : {np.max(vel_mag):.6f}")
            print(f"  Mean velocity: {np.mean(vel_mag):.6f}")
            print(f"  Max density  : {np.max(lb.rho):.6f}")

            plot_velocity_field(lb, step)

    end_time = time.time()

    # -------------------------------
    # Final reporting
    # -------------------------------
    if not converged:
        print("\nWARNING: Simulation did not fully converge")

    print("\n" + "=" * 60)
    print("Simulation completed")
    print("=" * 60)
    print(f"Final step        : {final_step}")
    print(f"Execution time    : {end_time - start_time:.2f} s")
    print(f"Convergence status: {'Yes' if converged else 'No'}")

    u_max = np.max(np.sqrt(lb.u[0]**2 + lb.u[1]**2))
    print(f"Maximum velocity : {u_max:.6f}")
    print(f"Reynolds number  : {Re:.2f}")

    # -------------------------------
    # FINAL POST-PROCESSING
    # -------------------------------
    print("\nPost-processing results...")

    # Velocity magnitude contours
    plot_velocity_field(lb, final_step)

    # Streamlines
    plot_streamlines_zoom(lb, final_step)

    # Pressure contours
    plot_pressure_field(lb, final_step)



if __name__ == "__main__":
    main()
