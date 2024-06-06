import numpy as np
import pandas as pd
import classes as cl

EPS = 1e-8

def find_collision_scatterer(pos, vel, scat_pos_r):
    rel_pos = pos - scat_pos_r[:2]
    a = np.dot(vel, vel)
    b = 2 * np.dot(vel, rel_pos)
    c = np.dot(rel_pos, rel_pos) - scat_pos_r[2]**2
    
    d = b**2 - 4*a*c

    # If d is negative, then particle trajectory doesn't intersect scatterer
    if d < 0:
        return None
    else:
        t1 = (-b + np.sqrt(d)) / (2*a)
        t2 = (-b - np.sqrt(d)) / (2*a)

        # Round small numbers to zero
        t1 = t1 if t1 > EPS else 0
        t2 = t2 if t2 > EPS else 0

        # Returns the smallest positive (non-zero) root
        t = min(t1, t2) if min(t1, t2) > 0 else max(t1, t2) if max(t1, t2) > 0 else None
        return t

def find_collisions_boundary(pos, vel, dimensions):
    w, h = dimensions
    max_x, min_x = w/2, -w/2
    max_y, min_y = h/2, -h/2

    t1 = (max_x - pos[0]) / vel[0] if vel[0] > 0 else None
    t2 = (min_x - pos[0]) / vel[0] if vel[0] < 0 else None
    t3 = (max_y - pos[1]) / vel[1] if vel[1] > 0 else None
    t4 = (min_y - pos[1]) / vel[1] if vel[1] < 0 else None

    t = [t1, t2, t3, t4]
    t = [x for x in t if x is not None and x >= 0] # >=0 since collision can happen on both scatterer and boundary

    return min(t) # Theoretically there should always be a collision

def check_particle_bounds(pos, dimensions):
    '''
    Checks if the particle is outside the boundary. If it is, it teleports the particle to the opposite side.
    Returns the new position of the particle and None if the particle is within the boundary.
    '''
    w, h = dimensions
    max_x, min_x = w/2, -w/2
    max_y, min_y = h/2, -h/2

    x = pos[0]
    y = pos[1]

    if x >= max_x:
        x = min_x
    elif x <= min_x:
        x = max_x
    if y >= max_y:
        y = min_y
    elif y <= min_y:
        y = max_y

    return np.array([x, y]) if x != pos[0] or y != pos[1] else None

def generate_tangential_trajectories(n, boundary, scatterers):
    '''
    Generates 2*s*n tangential trajectories for the system, where s is the number of scatterers.
    Each point generated has tangential trajectories in both directions.

    Note: even though the system has s scatterers, trajectories generated outside the boundary are not included.

    Parameters:
    n (int): Number of points to consider for each scatterer
    boundary (Boundary): Boundary object
    scatterers (list): List of Scatterer objects

    Returns:
    particles (list): List of Particle objects
    '''

    s = len(scatterers)
    particles = np.zeros((2*s*n, 4))
    thetas = np.linspace(0, 2*np.pi, n, endpoint=False)

    for i, scatterer in enumerate(scatterers):
        for j in range(n):
            x = scatterer.pos[0] + scatterer.radius*np.cos(thetas[j])
            y = scatterer.pos[1] + scatterer.radius*np.sin(thetas[j])
            vx = -np.sin(thetas[j])
            vy = np.cos(thetas[j])
            particles[i*n + j] = np.array([x, y, vx, vy])

    x_bounds = boundary.width/2
    y_bounds = boundary.height/2
    particles = particles[(particles[:, 0] >= -x_bounds) & (particles[:, 0] <= x_bounds) & (particles[:, 1] >= -y_bounds) & (particles[:, 1] <= y_bounds)]
    particles = particles[(particles[:, 2] != 0) | (particles[:, 3] != 0)] # Remove particles with zero velocity

    particles = [cl.Particle(arr=particles[i, :]) for i in range(len(particles))]
    return particles