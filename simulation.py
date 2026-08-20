import numpy as np
import config as cfg
from boundary_conditions import apply_x_bc, apply_y_bc, apply_pressure_bc


def run_simulation():
    cell_length = 1.0 / (cfg.N_POINTS_Y - 1)
    n_points_x = (cfg.N_POINTS_Y -1) * cfg.ASPECT_RATIO + 1
    x_range = np.linspace(0.0, 1.0 * cfg.ASPECT_RATIO, n_points_x)
    y_range = np.linspace(0.0, 1.0, cfg.N_POINTS_Y)

    coordinates_x, coordinates_y = np.meshgrid(x_range, y_range)

    #initial conditions

    velocity_x_prev = np.ones((cfg.N_POINTS_Y + 1, n_points_x))

    velocity_x_prev[:(cfg.BOX_HEIGHT_POINTS + 1), :] = 0

    #tom edge of the domain
    #(outside + inside) / 2 = 0  ===>  outside = -inside
    velocity_x_prev[-1, :] = - velocity_x_prev[-2, :]

    #top edge of box
    velocity_x_prev[cfg.BOX_HEIGHT_POINTS, 1:cfg.BOX_WITH_POINTS] = - velocity_x_prev[(cfg.BOX_HEIGHT_POINTS + 1), 1:cfg.BOX_WITH_POINTS]
    #right edge
    velocity_x_prev[1:(cfg.BOX_HEIGHT_POINTS + 1), cfg.BOX_WITH_POINTS] = 0.0

    #bottom of the domain
    #velocity_x_prev[0, :] = - velocity_x_prev[1, :]
    velocity_x_prev[0, (cfg.BOX_WITH_POINTS + 1):-1] = - velocity_x_prev[1, (cfg.BOX_WITH_POINTS + 1):-1]

    #values inside the box
    velocity_x_prev[:cfg.BOX_HEIGHT_POINTS, :cfg.BOX_WITH_POINTS] = 0.0

    velocity_y_prev = np.zeros((cfg.N_POINTS_Y, n_points_x + 1))

    pressure_prev = np.zeros((cfg.N_POINTS_Y + 1, n_points_x + 1))


    #prealloc arrays

    velocity_x_tent = np.zeros_like(velocity_x_prev)
    velocity_x_next = np.zeros_like(velocity_x_prev)

    velocity_y_tent = np.zeros_like(velocity_y_prev)
    velocity_y_next = np.zeros_like(velocity_y_prev)

    velocity_x_history, velocity_y_history = [], []

    for iter in range(cfg.N_TIME_STEPS):
        #update u
        diffusion_x = cfg.VISCOSITY * (
            (
                velocity_x_prev[1:-1, 2: ] + velocity_x_prev[2: , 1:-1] + velocity_x_prev[1:-1, :-2] + velocity_x_prev[:-2, 1:-1] - 4* velocity_x_prev[1:-1, 1:-1]
            ) / (cell_length**2)
        )
        convection_x = (
            (
                velocity_x_prev[1:-1, 2: ]**2 - velocity_x_prev[1: -1, :-2]**2
            ) / (2 * cell_length)
            +
            (
                velocity_y_prev[1: , 1:-2] + velocity_y_prev[1: , 2:-1] + velocity_y_prev[ :-1, 1:-2] + velocity_y_prev[ :-1, 2:-1]
            ) / 4
            *
            (
                velocity_x_prev[2: , 1:-1] - velocity_x_prev[ :-2, 1:-1]
            ) / (2 *cell_length)
        )
        pressure_grad_x = (
            (
                pressure_prev[1:-1, 2:-1] -  pressure_prev[1:-1, 1:-2]
            ) / (cell_length)
        )
    
        velocity_x_tent[1:-1, 1:-1] = (
            velocity_x_prev[1:-1, 1:-1] + cfg.TIME_STEP * (-pressure_grad_x + diffusion_x - convection_x)
        )
    
        velocity_x_tent = apply_x_bc(velocity_x_tent)
    
        #update v 
        diffusion_y = cfg.VISCOSITY * (
            (
                velocity_y_prev[1:-1, 2: ] + velocity_y_prev[2: , 1:-1] + velocity_y_prev[1:-1, :-2] + velocity_y_prev[:-2, 1:-1] - 4* velocity_y_prev[1:-1, 1:-1]
            ) / (cell_length**2)
            
        )
    
        convection_y = (
            (
                velocity_x_prev[2:-1, 1: ] + velocity_x_prev[2:-1, :-1] + velocity_x_prev[1:-2, 1: ] + velocity_x_prev[1:-2, :-1]
            ) / 4
            *
            (
                velocity_y_prev[1:-1, 2: ] - velocity_y_prev[1: -1, :-2]
            ) / (2 * cell_length)
            +
            (
                velocity_y_prev[2: , 1:-1]**2 - velocity_y_prev[ :-2, 1:-1]**2
            ) / (2 * cell_length)
        )
    
        pressure_grad_y = (
            (
                pressure_prev[2:-1, 1:-1] - pressure_prev[1:-2, 1:-1]
            ) / (cell_length)
        )
    
        velocity_y_tent[1:-1, 1:-1] = (
            velocity_y_prev[1:-1, 1:-1] + cfg.TIME_STEP * (- pressure_grad_y + diffusion_y - convection_y)
        )
    
        velocity_y_tent = apply_y_bc(velocity_y_tent)
    
        #divergance for pressure (at pressure nodes)
        divergence = (
            (
                velocity_x_tent[1:-1, 1: ] - velocity_x_tent[1:-1, :-1]
            ) / (cell_length)
            +
            (
                velocity_y_tent[1:, 1:-1] - velocity_y_tent[:-1, 1:-1]
            ) / (cell_length)
        )
    
        pressure_poisson_rhs = divergence / cfg.TIME_STEP
    
        #solve the pressure correction poission problem
    
        pressure_cor_prev = np.zeros_like(pressure_prev)
        for _ in range(cfg.N_PRESSURE_POISSION_ITERATIONS):
            pressure_cor_next = np.zeros_like(pressure_cor_prev)
            pressure_cor_next[1:-1, 1:-1] = 1/4 * ( 
                + pressure_cor_prev[1:-1, 2:] + pressure_cor_prev[2: , 1:-1] + pressure_cor_prev[1:-1, :-2] + pressure_cor_prev[:-2, 1:-1] 
                - cell_length**2 * pressure_poisson_rhs
                )
            #Apply BC: Homogeneous Neumann everywhere except for the left where homogeneous dirichlet
            pressure_cor_next = apply_pressure_bc(pressure_cor_next)
    
            pressure_cor_prev = pressure_cor_next
    
        pressure_next = pressure_prev + pressure_cor_next
    
    
        pressure_cor_grad_x = (
            (
                pressure_cor_next[1:-1, 2:-1] - pressure_cor_next[1:-1, 1:-2]
            ) / (cell_length)
        )
    
        velocity_x_next[1:-1, 1:-1] = (
            velocity_x_tent[1:-1, 1:-1] - cfg.TIME_STEP * pressure_cor_grad_x
        )
    
        pressure_cor_grad_y = (
            (
                pressure_cor_next[2:-1, 1:-1] - pressure_cor_next[1:-2, 1:-1]
            ) / (cell_length)
        )
    
        velocity_y_next[1:-1, 1:-1] = (
            velocity_y_tent[1:-1, 1:-1] - cfg.TIME_STEP * pressure_cor_grad_y
        )
    
        velocity_x_next = apply_x_bc(velocity_x_next)
        velocity_y_next = apply_y_bc(velocity_y_next)
    
    
        velocity_x_prev = velocity_x_next
        velocity_y_prev = velocity_y_next
        pressure_prev = pressure_next
    
        inflow_mass_rate_next = np.sum(velocity_x_next[1:-1, 0])
        outflow_mass_rate_next = np.sum(velocity_x_next[1:-1, -1])
    
        # Visualization
        if iter % cfg.PLOT_EVERY == 0:
            print(f"Inflow: {inflow_mass_rate_next}, Outflow: {outflow_mass_rate_next}")
            velocity_x_vertex_centered = (
                (
                    velocity_x_next[1: , :] + velocity_x_next[:-1, :]
                ) / 2
            )
            velocity_y_vertex_centered = (
                (
                    velocity_y_next[:, 1: ] + velocity_y_next[:, :-1]
                ) / 2
            )
            
            velocity_x_vertex_centered[:(cfg.BOX_HEIGHT_POINTS + 1), :(cfg.BOX_WITH_POINTS + 1)] = 0.0
            velocity_y_vertex_centered[:(cfg.BOX_HEIGHT_POINTS + 1), :(cfg.BOX_WITH_POINTS + 1)] = 0.0
    
            velocity_x_history.append(velocity_x_vertex_centered)
            velocity_y_history.append(velocity_y_vertex_centered)

    return np.array(velocity_x_history), np.array(velocity_y_history), coordinates_x, coordinates_y