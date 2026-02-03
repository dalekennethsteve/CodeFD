"""
main_rotating_simple.py - Fixed version with single loop
"""
import numpy as np
import time
from rotating_cylinder import RotatingCylinder
from visualization_rotating import plot_streamlines_rotating, plot_streamlines_zoom_rotating, \
    plot_streamlines_with_velocity_color, plot_velocity_rotating, plot_pressure_rotating, plot_force_history
from parameters import *
from utils import check_convergence


def main():
    print("Rotating Cylinder Simulation with Force Calculation")
    print("=" * 50)

    # Set rotation speed
    omega_rot = 0.065  # Change this value for different rotations

    print(f"Angular velocity: ω = {omega_rot}")
    print(f"Cylinder at ({cylinder_x}, {cylinder_y}), radius {cylinder_radius}")

    # Initialize simulation
    lb = RotatingCylinder(omega_rot=omega_rot)

    # Simulation parameters
    total_steps = 10000
    plot_interval = 500
    convergence_check_interval = 100
    convergence_tolerance = 1e-6

    # FORCE TRACKING
    force_history = []  # Store (step, drag, lift) tuples
    force_check_interval = 50  # Calculate forces every N steps

    print(f"\nRunning {total_steps} steps...")
    start_time = time.time()

    # Variables for tracking
    converged = False
    max_change = 0
    final_step = total_steps

    # SINGLE LOOP - everything happens here
    for step in range(total_steps):
        # Store velocity every 100 steps for convergence check
        if step % convergence_check_interval == 0:
            prev_u = lb.u.copy()

        # Run LBM step
        lb.step()

        # Calculate forces at intervals
        if step % force_check_interval == 0 and step > 0:
            drag_force, lift_force = lb.calculate_forces()
            force_history.append((step, drag_force, lift_force))
            if len(force_history) % 5 == 0:  # Print every 5 force calculations
                print(f"Step {step}: Drag = {drag_force:.6f}, Lift = {lift_force:.6f}")

        # Check convergence every 100 steps after step 0
        if step % convergence_check_interval == 0 and step > 0:
            converged, max_change = check_convergence(lb.u, prev_u, convergence_tolerance)
            if converged:
                final_step = step
                print(f"\n✓ CONVERGENCE ACHIEVED at step {step}")
                print(f"  Maximum velocity change: {max_change:.2e}")
                break

        # Plot progress at intervals
        if step % plot_interval == 0:
            print(f"Step {step}/{total_steps} (max change: {max_change:.2e})")
            #plot_velocity_rotating(lb, step, omega_rot)
            #plot_streamlines_rotating(lb, step, omega_rot)
            plot_streamlines_zoom_rotating(lb, step, omega_rot)

    end_time = time.time()

    # Handle non-convergence
    if not converged:
        print(f"\n⚠ WARNING: Did not converge within {total_steps} steps")
        print(f"  Final maximum velocity change: {max_change:.2e}")
        final_step = total_steps

    # Calculate final forces
    final_drag, final_lift = lb.calculate_forces()

    # Final plots
    print(f"\nSimulation completed in {end_time - start_time:.2f} seconds")
    print(f"Final step: {final_step}")

    plot_velocity_rotating(lb, final_step, omega_rot)
    plot_pressure_rotating(lb, final_step, omega_rot)


    # NEW: Plot force history
    plot_force_history(force_history, omega_rot)
    plot_streamlines_zoom_rotating(lb, step, omega_rot=0.065, save_path='output')

    # Show final statistics
    vel_mag = lb.get_velocity_magnitude()
    max_vel = np.max(vel_mag)
    avg_vel = np.mean(vel_mag[~lb.obstacle])

    print(f"\n=== FINAL STATISTICS ===")
    print(f"Velocity:")
    print(f"  Max velocity: {max_vel:.6f}")
    print(f"  Average velocity in fluid: {avg_vel:.6f}")
    print(f"  Tangential velocity at surface: {abs(omega_rot * cylinder_radius):.4f}")

    print(f"\nForces on cylinder:")
    print(f"  Drag force: {final_drag:.6f} (lattice units)")
    print(f"  Lift force: {final_lift:.6f} (lattice units)")
    print(f"  Lift/Drag ratio: {abs(final_lift / final_drag):.3f}" if final_drag != 0 else "  Lift/Drag ratio: ∞")

    print(f"\nPerformance:")
    print(f"  Total steps: {final_step}")
    print(f"  Convergence: {'YES' if converged else 'NO'}")
    print(f"  Final max change: {max_change:.2e}")

    # Convert to dimensionless coefficients
    U_ref = abs(omega_rot * cylinder_radius)  # Reference velocity = surface speed
    if U_ref > 0:
        # Simple approximation for coefficients
        C_d = 2 * final_drag / (U_ref ** 2 * cylinder_radius)  # Drag coefficient
        C_l = 2 * final_lift / (U_ref ** 2 * cylinder_radius)  # Lift coefficient
        print(f"\nDimensionless coefficients:")
        print(f"  Drag coefficient C_d ≈ {C_d:.4f}")
        print(f"  Lift coefficient C_l ≈ {C_l:.4f}")

if __name__ == "__main__":
    main()