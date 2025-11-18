import numpy as np
import time
from lattice_boltzmann import LatticeBoltzmann
from boundary_conditions import apply_boundary_conditions
from parameters import tau, nu, body_force
from src.utils import calculate_Reynoldsnumber
from visualization import plot_velocity_field
from utils import check_convergence


def main():
    # Initialize LBM simulation
    print("Initializing LBM simulation with periodic boundary conditions...")
    lb = LatticeBoltzmann()

    # Simulation parameters
    total_steps = 10000
    plot_interval = 500
    convergence_check_interval = 100
    convergence_tolerance = 1e-6
    Re = calculate_Reynoldsnumber(lb)

    print(f"Running simulation for {total_steps} steps...")


    # Check initial conditions
    print(f"Initial max velocity: {np.max(np.abs(lb.u)):.6f}")
    print(f"Initial max density: {np.max(lb.rho):.6f}")

    start_time = time.time()
    converged = False
    final_step = total_steps

    for step in range(total_steps):
        # Store previous velocity for convergence check
        if step % convergence_check_interval == 0:
            prev_u = lb.u.copy()

        # Main LBM steps
        lb.step()

        # Check convergence
        if step % convergence_check_interval == 0 and step > 0:
            converged, max_change = check_convergence(lb.u, prev_u, convergence_tolerance)

            if converged:
                final_step = step
                print(f"Convergence achieved at step {step}")
                print(f"Maximum velocity change: {max_change:.2e}")
                break

        # Debug output
        if step % plot_interval == 0:
            max_vel = np.max(np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2))
            avg_vel = np.mean(np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2))

            print(f"Step {step}/{total_steps}")
            print(f"  Max velocity: {max_vel:.6f}")
            print(f"  Avg velocity: {avg_vel:.6f}")
            print(f"  Max density: {np.max(lb.rho):.6f}")

            plot_velocity_field(lb, step)

    if not converged:
        print(f"Warning: Simulation did not converge within {total_steps} steps")
        print(f"Final maximum velocity change: {max_change:.2e}")

    end_time = time.time()
    print(f"Simulation completed in {end_time - start_time:.2f} seconds")
    print(f"Final step: {final_step}")

    # Final plots
    plot_velocity_field(lb, final_step)

    # Print final statistics
    u_max = np.max(np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2))
    print(f"\nFinal statistics:")
    print(f"Maximum velocity: {u_max:.6f}")
    print(f"Reynolds number: {Re:.2f}")
    print(f"Convergence status: {'Yes' if converged else 'No'}")


if __name__ == "__main__":
    main()