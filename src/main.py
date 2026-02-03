import time
import numpy as np
from lattice_boltzmann import LatticeBoltzmann
from visualization import plot_streamlines_zoom
from utils import check_convergence
from parameters import *
import matplotlib.pyplot as plt


def main():
    # Initialize LBM simulation
    print("Initializing LBM simulation...")
    lb = LatticeBoltzmann()

    # Count obstacle nodes
    obstacle_nodes = np.sum(lb.obstacle)
    print(f"Obstacle nodes: {obstacle_nodes}")

    # Simulation parameters
    total_steps = 10000
    plot_interval = 1000
    convergence_check_interval = 100
    convergence_tolerance = 1e-6

    print(f"Running simulation for {total_steps} steps...")

    # Check initial conditions
    print(f"Initial max velocity: {np.max(np.abs(lb.u)):.6f}")
    print(f"Initial max density: {np.max(lb.rho):.6f}")

    start_time = time.time()
    converged = False
    final_step = total_steps
    max_change = 0

    # Track whether we should save final plots
    save_final_plots = True

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

        # Debug output and plotting - ONLY VELOCITY AND PRESSURE EVERY 1000 STEPS
        if step % plot_interval == 0:
            max_vel = np.max(np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2))
            avg_vel = np.mean(np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2))

            print(f"Step {step}/{total_steps}")
            print(f"  Max velocity: {max_vel:.6f}")
            print(f"  Avg velocity: {avg_vel:.6f}")
            print(f"  Max density: {np.max(lb.rho):.6f}")

            #plot_velocity_field(lb, step)
            #plot_pressure_field(lb, step)
            #plot_streamlines_with_velocity(lb, step)
            plot_streamlines_zoom(lb, step)

            # Clear figures to prevent accumulation
            plt.close('all')

    if not converged:
        print(f"Warning: Simulation did not converge within {total_steps} steps")
        print(f"Final maximum velocity change: {max_change:.2e}")

    end_time = time.time()
    print(f"Simulation completed in {end_time - start_time:.2f} seconds")
    print(f"Final step: {final_step}")

    # Print final statistics
    u_max = np.max(np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2))
    print(f"\nFinal statistics:")
    print(f"Maximum velocity: {u_max:.6f}")
    print(f"Convergence status: {'Yes' if converged else 'No'}")
    print(f"Reynolds number: {Re:0.2f}")

    # Plot final results - ALL THREE PLOTS AT THE END
    if save_final_plots:
        print("\nGenerating final plots...")

        # Create a figure with subplots for final results
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # Final velocity plot
        print("  Generating final velocity plot...")
        X, Y = np.meshgrid(np.arange(lb.nx), np.arange(lb.ny), indexing='ij')
        velocity_magnitude = np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2)

        im1 = axes[0].contourf(X, Y, velocity_magnitude, levels=50, cmap='jet')
        plt.colorbar(im1, ax=axes[0], label='Velocity Magnitude')
        axes[0].set_title(f'Final Velocity - Step {final_step}')
        axes[0].set_xlabel('X')
        axes[0].set_ylabel('Y')
        axes[0].axis('equal')


        # Add square obstacle to all plots
        square_half_size = cylinder_radius
        square_x = cylinder_x - square_half_size
        square_y = cylinder_y - square_half_size
        square_width = 2 * square_half_size

        from matplotlib.patches import Rectangle
        square = Rectangle((square_x, square_y), square_width, square_width,
                           linewidth=2, edgecolor='black', facecolor='none')
        axes[0].add_patch(square)


        ''''# Final pressure plot
        print("  Generating final pressure plot...")
        pressure = (1 / 3) * lb.rho
        im2 = axes[1].contourf(X, Y, pressure, levels=50, cmap='jet')
        plt.colorbar(im2, ax=axes[1], label='Pressure')
        axes[1].set_title(f'Final Pressure - Step {final_step}')
        axes[1].set_xlabel('X')
        axes[1].set_ylabel('Y')
        axes[1].axis('equal')
        
        axes[1].add_patch(Rectangle((square_x, square_y), square_width, square_width,
                                    linewidth=2, edgecolor='black', facecolor='none'))
      

        # Final streamlines plot
        print("  Generating final streamlines plot...")
        x = np.arange(0, lb.nx, 1)
        y = np.arange(0, lb.ny, 1)
        axes[2].streamplot(x, y, lb.u[0].T, lb.u[1].T,
                           color='darkblue',
                           linewidth=0.6,
                           density=2,
                           arrowsize=0.6 )
        axes[2].set_title(f'Final Streamlines - Step {final_step}')
        axes[2].set_xlabel('X')
        axes[2].set_ylabel('Y')
        axes[2].axis('equal')
        
        axes[2].add_patch(Rectangle((square_x, square_y), square_width, square_width,
                                    linewidth=2, edgecolor='black', facecolor='black', alpha=0.8))
        

        plt.tight_layout()
        plt.show()'''

        # Also show individual full-size plots
        print("\nGenerating individual full-size plots...")
        #plot_velocity_field(lb, final_step)
        #plot_pressure_field(lb, final_step)
        #plot_streamlines_with_velocity(lb, final_step)
        #plot_velocity_with_streamlines(lb, final_step)
        plot_streamlines_zoom(lb, final_step)

    # Keep plots open at the end
    plt.show()


if __name__ == "__main__":
    main()