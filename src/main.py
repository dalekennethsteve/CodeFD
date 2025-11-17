import numpy as np
import time
from lattice_boltzmann import LatticeBoltzmann
from boundary_conditions import apply_boundary_conditions
from src.parameters import tau, nu
from visualization import plot_velocity_field



def main():
    # Initialize LBM simulation
    print("Initializing LBM simulation...")
    lb = LatticeBoltzmann()

    # Simulation parameters
    total_steps = 1000
    plot_interval = 200
    convergence_check_interval = 100

    print(f"Running simulation for {total_steps} steps...")
    print(f"Domain size: {lb.nx} x {lb.ny}")
    print(f"Relaxation time tau: {tau}")
    print(f"Kinematic viscosity: {nu:.6f}")

    # Check initial conditions
    print(f"Initial max velocity: {np.max(np.abs(lb.u)):.6f}")
    print(f"Initial max density: {np.max(lb.rho):.6f}")

    start_time = time.time()

    for step in range(total_steps):
        # Store previous velocity for convergence check
        if step % convergence_check_interval == 0:
            prev_u = lb.u.copy()

        # Main LBM steps
        lb.stream()
        apply_boundary_conditions(lb)
        lb.collide()

        # Debug output
        if step % plot_interval == 0:
            max_vel = np.max(np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2))
            avg_vel = np.mean(np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2))

            print(f"Step {step}/{total_steps}")
            print(f"  Max velocity: {max_vel:.6f}")
            print(f"  Avg velocity: {avg_vel:.6f}")
            print(f"  Max density: {np.max(lb.rho):.6f}")


            plot_velocity_field(lb, step)


    end_time = time.time()
    print(f"Simulation completed in {end_time - start_time:.2f} seconds")

    # Final plots
    plot_velocity_field(lb, total_steps)


    # Print final statistics
    u_max = np.max(np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2))
    print(f"\nFinal statistics:")
    print(f"Maximum velocity: {u_max:.6f}")



if __name__ == "__main__":
    main()