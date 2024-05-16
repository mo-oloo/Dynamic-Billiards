import numpy as np
import src.classes as cl

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

