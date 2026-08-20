# Fluid Flow Simulation

A 2D computational fluid dynamics (CFD) simulation written in Python. This project solves the Navier-Stokes equations using the finite difference method obtained from taylor's expensions to simulate fluid flow around a rectangular obstacle. 

The simulation outputs animated GIFs visualizing the flow using both contour maps and stream plots.

## Project Structure

* **`config.py`**: Contains all tunable simulation parameters (grid size, viscosity, time steps, etc.).
* **`boundary_conditions.py`**: Handles the boundary logic for velocity and pressure arrays.
* **`simulation.py`**: The core physics engine. It runs the main time-stepping loop and calculates pressure corrections and momentum.
* **`main.py`**: The entry point of the simulation. It runs the simulation and generates the Matplotlib animations.

## Requirements

You will need Python 3 installed, along with a few standard python libraries like numpy, matplotlib etc.. 

```bash
pip install numpy matplotlib
