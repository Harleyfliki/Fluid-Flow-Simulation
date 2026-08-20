import numpy as np
from config import BOX_HEIGHT_POINTS, BOX_WITH_POINTS

def apply_x_bc(v):
    #inflow
    #v[1:-1, 0] = 1.0
    v[(BOX_HEIGHT_POINTS + 1):-1, 0] = 1

    #outflow
    inflow_mass_rate = np.sum(v[(BOX_HEIGHT_POINTS + 1):-1, 0])
    outflow_mass_rate = np.sum(v[1:-1, -2])
    #v[1:-1, -1] = v[1:-1, -2]
    v[1:-1, -1] = v[1:-1, -2] * inflow_mass_rate / outflow_mass_rate

    #top of box
    v[BOX_HEIGHT_POINTS, 1:BOX_WITH_POINTS] = - v[BOX_HEIGHT_POINTS + 1, 1:BOX_WITH_POINTS]

    #right edge of the box
    v[1:(BOX_HEIGHT_POINTS + 1), BOX_WITH_POINTS] = 0.0

    #bottom edge of domain
    v[0, (BOX_WITH_POINTS +  1):-1] = - v[1, (BOX_WITH_POINTS + 1):-1]

    #top edge of domain
    #v[0, :] = - v[1, :]
    v[-1, :] = -v[-2, :]

    #set all U= 0 inside the box
    v[:BOX_HEIGHT_POINTS, :BOX_WITH_POINTS] = 0

    return v

def apply_y_bc(v): 
    #inflow
    #v[1:-1, 0] = - v[1:-1, 1]
    v[(BOX_HEIGHT_POINTS + 1):-1, 0] = -v[(BOX_HEIGHT_POINTS + 1):-1, 1]

    #outflow
    # (outside - inside)  * dx = 0 Neumann BC    outside = inside
    v[1:-1, -1] = v[1:-1, -2]

    #Top edge of box
    v[BOX_HEIGHT_POINTS, 1:(BOX_WITH_POINTS + 1)] = 0.0

    #right edge of the box
    v[1:(BOX_HEIGHT_POINTS + 1), BOX_WITH_POINTS] = -v[1:(BOX_HEIGHT_POINTS + 1), (BOX_WITH_POINTS + 1)]


    #bottom edge of domain
    v[0, (BOX_WITH_POINTS + 1):] = 0.0
    # v[0, :] = 0.0

    #top edge of domain
    v[-1, :] = 0.0

    #set all v= 0 inside the box
    v[:BOX_HEIGHT_POINTS, :BOX_WITH_POINTS] = 0

    return v

def apply_pressure_bc(p):
    #inflow
    # p[1:-1, 0] = p[1:-1, 1]
    p[(BOX_HEIGHT_POINTS + 1):-1, 0] = p[(BOX_HEIGHT_POINTS + 1):-1, 1]

    #outflow
    p[1:-1, -1] = - p[1:-1, -2]

    #top edge of box
    p[BOX_HEIGHT_POINTS, 1:(BOX_WITH_POINTS + 1)] = p[(BOX_HEIGHT_POINTS + 1), 1:(BOX_WITH_POINTS + 1)]

    # right edge of box
    p[1:(BOX_HEIGHT_POINTS + 1), BOX_WITH_POINTS] = p[1:(BOX_HEIGHT_POINTS + 1), (BOX_WITH_POINTS + 1)]

    #bottom edge of domain
    p[0, (BOX_WITH_POINTS + 1):-1] = p[1, (BOX_WITH_POINTS + 1):-1]
    # p[0, :] = p[1, :]

    #top edge of domain
    p[-1, :] = p[-2, :]

    #set all cor values inside box to 0
    p[:BOX_HEIGHT_POINTS, :BOX_WITH_POINTS] = 0.0

    return p