import numpy as np
import matplotlib.subplots as plt
from matplotlib.animation import FuncAnimation
from simulation import run_simulation

def main():
    print("Starting fluid flow simulation...")

    velocity_x_history, velocity_y_history, coordinates_x, coordinates_y = run_simulation()

    print("Simulation complete. Generating animations...")
  
    plt.style.use('dark_background') 
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_title('Fluid Flow Simulation (Stream Plot)')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    # Initial stream plot
    stream = ax.streamplot(
        coordinates_x, coordinates_y,
        velocity_x_history[0], velocity_y_history[0],
        color=np.sqrt(velocity_x_history[0]**2 + velocity_y_history[0]**2),
        cmap='viridis', linewidth=2, density=2
    )
    fig.colorbar(
        plt.cm.ScalarMappable(cmap='viridis'),
        ax=ax,
        label='Velocity Magnitude'
    )

    def update(frame):
        """Update function for stream plot animation."""
        ax.clear()
        ax.set_title('Fluid Flow Simulation (Stream Plot)')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        
        stream = ax.streamplot(
            coordinates_x, coordinates_y,
            velocity_x_history[frame], velocity_y_history[frame],
            color=np.sqrt(velocity_x_history[frame]**2 + velocity_y_history[frame]**2),
            cmap='viridis', linewidth=2, density=2
        )
        return stream.lines

    # Create the animation
    anim = FuncAnimation(fig, update, frames=len(velocity_x_history), interval=200, blit=False)

    anim.save("fluid_flow_streamplot.gif", fps=5)
    plt.show()

if __name__ == "__main__":
    main()