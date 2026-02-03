import time
import numpy as np
from lattice_boltzmann import LatticeBoltzmann
from utils import check_convergence
from parameters import *
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


def plot_beautiful_streamlines(lb, step, show_profile=True):
    """Plot beautiful streamlines colored by velocity magnitude"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={'width_ratios': [2, 1]})

    # Create grid for streamplot
    X, Y = np.meshgrid(np.arange(lb.nx), np.arange(lb.ny), indexing='ij')

    # Calculate velocity magnitude for coloring
    velocity_magnitude = np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2)

    # Create colormap
    cmap = plt.cm.jet
    norm = Normalize(vmin=0, vmax=np.max(velocity_magnitude) * 1.1)

    # Subsample for cleaner streamlines (every 2 points)
    subsample = 2
    X_sub = X[::subsample, ::subsample]
    Y_sub = Y[::subsample, ::subsample]
    U_sub = lb.u[0, ::subsample, ::subsample]
    V_sub = lb.u[1, ::subsample, ::subsample]
    Vmag_sub = velocity_magnitude[::subsample, ::subsample]

    # Plot beautiful streamlines with velocity magnitude coloring
    stream = ax1.streamplot(X_sub.T, Y_sub.T, U_sub.T, V_sub.T,
                            color=Vmag_sub.T,
                            cmap=cmap,
                            norm=norm,
                            linewidth=1.5,
                            density=2.5,
                            arrowsize=1.2,
                            arrowstyle='->',
                            minlength=0.1,
                            maxlength=4.0)

    # Add walls
    ax1.axhline(y=0, color='black', linewidth=3, alpha=0.9)
    ax1.axhline(y=ny - 1, color='black', linewidth=3, alpha=0.9)

    # Add colorbar
    cbar = plt.colorbar(stream.lines, ax=ax1, orientation='vertical', pad=0.02)
    cbar.set_label('Velocity Magnitude', fontsize=12)

    # Beautify the plot
    ax1.set_title(f'Poiseuille Flow - Step {step}', fontsize=14, fontweight='bold', pad=15)
    ax1.set_xlabel('X (Flow Direction)', fontsize=12)
    ax1.set_ylabel('Y (Channel Height)', fontsize=12)
    ax1.set_xlim(0, nx)
    ax1.set_ylim(0, ny)
    ax1.grid(True, alpha=0.2, linestyle='--')
    ax1.set_aspect('equal', adjustable='box')

    # Add flow information
    max_vel = np.max(velocity_magnitude)
    avg_vel = np.mean(np.abs(lb.u[0]))
    error = abs(max_vel - u_max_predicted) / u_max_predicted * 100

    info_text = f'Max velocity: {max_vel:.4f}\nPredicted: {u_max_predicted:.4f}\nError: {error:.2f}%\nRe: {Re:.1f}'

    ax1.text(0.02, 0.98, info_text,
             transform=ax1.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))

    # Right subplot: Velocity profile
    if show_profile:
        # Get velocity profile at middle of channel
        center_x = nx // 2
        velocity_profile = lb.u[0, center_x, :]
        y_positions = np.arange(ny)

        # Plot velocity profile with gradient color
        for i in range(ny - 1):
            y_vals = [y_positions[i], y_positions[i + 1]]
            u_vals = [velocity_profile[i], velocity_profile[i + 1]]
            ax2.plot(u_vals, y_vals, color=cmap(norm(np.mean(u_vals))), linewidth=3)

        # Theoretical parabolic profile (dashed line)
        theoretical_profile = u_max_predicted * (1 - ((y_positions - ny / 2) / (ny / 2)) ** 2)
        ax2.plot(theoretical_profile, y_positions, 'k--', linewidth=2,
                 label='Theoretical', alpha=0.7)

        ax2.set_xlabel('Velocity', fontsize=12)
        ax2.set_ylabel('Y Position', fontsize=12)
        ax2.set_title('Velocity Profile at Center', fontsize=14, fontweight='bold', pad=15)
        ax2.legend(loc='upper right', fontsize=11)
        ax2.grid(True, alpha=0.3, linestyle='--')
        ax2.set_xlim(0, max_vel * 1.2)
        ax2.set_ylim(0, ny)

        # Add velocity magnitude colorbar for profile
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar2 = plt.colorbar(sm, ax=ax2, orientation='vertical', pad=0.02)
        cbar2.set_label('Velocity Magnitude', fontsize=12)

    plt.tight_layout()
    plt.show()


def main():
    # Initialize LBM simulation
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║               POISEUILLE FLOW SIMULATION                    ║")
    print("╚══════════════════════════════════════════════════════════════╝")
    print(f"\n📐 Channel dimensions: {nx} × {ny}")
    print(f"⚙️  Relaxation time (τ): {tau}")
    print(f"🌊 Kinematic viscosity (ν): {nu:.6f}")
    print(f"📈 Body force (F_x): {body_force[0]:.6f}")
    print(f"🎯 Predicted max velocity: {u_max_predicted:.6f}")
    print(f"🌀 Reynolds number: {Re:.2f}")
    print("-" * 60)

    lb = LatticeBoltzmann()

    # Simulation parameters
    total_steps = 10000
    plot_interval = 1000
    convergence_check_interval = 200
    convergence_tolerance = 1e-8

    print(f"\n🚀 Running simulation for {total_steps} steps...")
    print(f"📊 Plotting every {plot_interval} steps")
    print("-" * 60)

    start_time = time.time()
    converged = False
    final_step = total_steps
    max_change = 0
    velocity_history = []

    for step in range(total_steps):
        # Store previous velocity for convergence check
        if step % convergence_check_interval == 0:
            prev_u = lb.u.copy()

        # Main LBM steps
        lb.step()

        # Track maximum velocity for monitoring
        if step % 100 == 0:
            max_vel = np.max(np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2))
            velocity_history.append((step, max_vel))

        # Check convergence
        if step % convergence_check_interval == 0 and step > 0:
            converged, max_change = check_convergence(lb.u, prev_u, convergence_tolerance)

            if converged:
                final_step = step
                print(f"\n✅ Convergence achieved at step {step}")
                print(f"   Maximum velocity change: {max_change:.2e}")
                break

        # Plot beautiful streamlines at intervals
        if step % plot_interval == 0:
            max_vel = np.max(np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2))
            avg_vel = np.mean(np.abs(lb.u[0]))
            error = abs(max_vel - u_max_predicted) / u_max_predicted * 100

            print(f"\n📈 Step {step:6d}/{total_steps}")
            print(f"   Max velocity: {max_vel:.6f} (expected: {u_max_predicted:.6f})")
            print(f"   Average velocity: {avg_vel:.6f}")
            print(f"   Relative error: {error:.2f}%")

            # Plot beautiful streamlines with velocity profile
            plot_beautiful_streamlines(lb, step, show_profile=(step > 0))

            # Clear figure to prevent accumulation
            plt.close('all')

    if not converged:
        print(f"\n⚠️  Warning: Simulation did not converge within {total_steps} steps")
        print(f"   Final maximum velocity change: {max_change:.2e}")

    end_time = time.time()
    print("\n" + "=" * 60)
    print("🎉 SIMULATION COMPLETE")
    print("=" * 60)
    print(f"⏱️  Total time: {end_time - start_time:.2f} seconds")
    print(f"🔢 Final step: {final_step}")
    print(f"⚡ Steps per second: {final_step / (end_time - start_time):.0f}")

    # Final statistics
    u_max = np.max(np.sqrt(lb.u[0] ** 2 + lb.u[1] ** 2))
    u_center = lb.u[0, nx // 2, ny // 2]
    error = abs(u_max - u_max_predicted) / u_max_predicted * 100

    print(f"\n📊 FINAL RESULTS:")
    print(f"   Max velocity (simulated):  {u_max:.6f}")
    print(f"   Max velocity (theoretical): {u_max_predicted:.6f}")
    print(f"   Centerline velocity:       {u_center:.6f}")
    print(f"   Relative error:            {error:.2f}%")
    print(f"   Convergence:               {'✅ CONVERGED' if converged else '❌ NOT CONVERGED'}")
    print(f"   Reynolds number:           {Re:.2f}")

    # Plot final beautiful results
    print("\n🎨 Generating final beautiful plot...")

    # Final detailed streamline plot
    plot_beautiful_streamlines(lb, final_step, show_profile=True)

    # Additional plot: Velocity evolution over time
    if velocity_history:
        steps, velocities = zip(*velocity_history)

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(steps, velocities, 'b-', linewidth=2.5, label='Simulated')
        ax.axhline(y=u_max_predicted, color='r', linestyle='--',
                   linewidth=2, label='Theoretical', alpha=0.7)

        # Fill between for visual appeal
        ax.fill_between(steps, 0, velocities, alpha=0.3, color='blue')

        # Add convergence line if converged
        if converged:
            ax.axvline(x=final_step, color='g', linestyle=':',
                       linewidth=2, alpha=0.7, label='Convergence')

        ax.set_xlabel('Time Step', fontsize=12)
        ax.set_ylabel('Maximum Velocity', fontsize=12)
        ax.set_title('Velocity Evolution During Simulation', fontsize=14, fontweight='bold')
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_ylim(0, max(velocities) * 1.1)

        # Add statistics box
        stats_text = f'Final velocity: {u_max:.4f}\nError: {error:.2f}%\nSteps to converge: {final_step}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                verticalalignment='top', fontsize=11,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

        plt.tight_layout()
        plt.show()

    print("\n✨ Simulation complete. Close plots to exit.")
    print("=" * 60)

    # Keep plots open
    plt.show()


if __name__ == "__main__":
    main()