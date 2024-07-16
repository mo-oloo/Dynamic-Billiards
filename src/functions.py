import numpy as np
import pandas as pd
import classes as cl

EPS = 1e-6

def find_collision_scatterer(pos, vel, scat_pos_r):
    rel_pos = pos - scat_pos_r[:2]
    a = np.dot(vel, vel)
    b = 2 * np.dot(vel, rel_pos)
    c = np.dot(rel_pos, rel_pos) - scat_pos_r[2]**2
    
    d = b**2 - 4*a*c

    # If d is negative, then particle trajectory doesn't intersect scatterer
    if d < 0 or a == 0:
        return None
    else:
        d = np.sqrt(d)
        t1 = (-b + d) / (2*a)
        t2 = (-b - d) / (2*a)

        # Round small numbers to zero
        t1 = t1 if t1 > EPS else 0
        t2 = t2 if t2 > EPS else 0

        # Returns the smallest positive (non-zero) root
        positive_roots = [t for t in [t1, t2] if t > 0]
        t = min(positive_roots) if positive_roots else None
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
    change_vector = np.array([0, 0])

    if x >= max_x:
        x = min_x
        change_vector[0] += 1
    elif x <= min_x:
        x = max_x
        change_vector[0] -= 1
    if y >= max_y:
        y = min_y
        change_vector[1] += 1
    elif y <= min_y:
        y = max_y
        change_vector[1] -= 1

    if x != pos[0] or y != pos[1]:
        return np.array([x, y]), change_vector
    else:
        return None, None

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
    particles = []
    scat = []
    thet = []
    incvecs = []
    thetas = np.linspace(0, 2*np.pi, n, endpoint=True)

    for i, scatterer in enumerate(scatterers):
        for j in range(n):
            x = scatterer.pos[0] + scatterer.radius*np.cos(thetas[j])
            y = scatterer.pos[1] + scatterer.radius*np.sin(thetas[j])
            vx = -np.sin(thetas[j])
            vy = np.cos(thetas[j])
            # particles[i*n + j] = np.array([x, y, vx, vy])
            # particles[i*n + j + n] = np.array([x, y, -vx, -vy])
            particles.extend((np.array([x, y, vx, vy]), np.array([x, y, -vx, -vy])))
            scat.extend((i, i))
            thet.extend((thetas[j], thetas[j]))
            incvecs.extend((np.pi/2, -np.pi/2))
    particles = np.array(particles)
    scat = np.array(scat)
    thet = np.array(thet)
    incvecs = np.array(incvecs)
    
    x_bounds = boundary.width/2
    y_bounds = boundary.height/2

    # Get index of particles within the boundary
    in_bounds = (particles[:, 0] >= -x_bounds) & (particles[:, 0] <= x_bounds) & (particles[:, 1] >= -y_bounds) & (particles[:, 1] <= y_bounds)
    particles = particles[in_bounds]
    scat = scat[in_bounds]
    thet = thet[in_bounds]
    incvecs = incvecs[in_bounds]

    particles = [cl.Particle(arr=particles[i, :]) for i in range(len(particles))]
    return particles, scat, thet, incvecs

def generate_random_trajectories(size, scatterer, arange=False):
    '''
    Generates a grid of initial conditions for a particle on the center scatterer.

    Parameters:
    size (array-like): Number of points to generate for the theta and phi values or the spacing between the points if arange=True
    scatterer (Scatterer): Scatterer object
    arange (bool): Whether to use np.arange or np.linspace. If True, size is the spacing between the points.

    Returns:
    particles (list): List of Particle objects
    '''
    m, n = size
    if arange:
        thetas = np.arange(0, 2*np.pi, m)
        phis = np.arange(-np.pi/2, np.pi/2, n)
    else:
        thetas = np.linspace(0, 2*np.pi, m, endpoint=True)
        phis = np.linspace(-np.pi/2, np.pi/2, n, endpoint=True)

    particles = []
    scat = []
    thet = []
    vel = []

    x0 = scatterer.pos[0]
    y0 = scatterer.pos[1]

    for theta in thetas:
        for phi in phis:
            x = x0 + scatterer.radius*np.cos(theta)
            y = y0 + scatterer.radius*np.sin(theta)
            vx = np.cos(theta + phi)
            vy = np.sin(theta + phi)
            particles.append(np.array([x, y, vx, vy]))
            scat.append(0)
            thet.append(theta)
            vel.append(phi)
    particles = np.array(particles)
    scat = np.array(scat)
    thet = np.array(thet)
    vel = np.array(vel)

    particles = [cl.Particle(arr=particles[i, :]) for i in range(len(particles))]

    return particles, scat, thet, vel

def process_trajectory_sequence(data, big_scat, a, b):
    data = data.dropna(subset=['symbol'])
    data = data.reset_index(drop=True)
    # Convert each symbol to a tuple
    data.loc[:, 'symbol'] = data['symbol'].apply(lambda x: tuple(x))
    return check_over_under(data, big_scat, a, b)
    # return data

def check_over_under(data, big_scat, a, b):
    for i in range(1, len(data)):
        if data['scatterer_hit'].iloc[i-1] == data['scatterer_hit'].iloc[i] and data['scatterer_hit'].iloc[i] == big_scat:
            symbol = data['symbol'].iloc[i]
            symbol_abs = tuple(abs(x) for x in symbol) # Get absolute value
            if symbol_abs != (2*a, 2*b):
                continue
            data2 = data.iloc[i]
            data1 = data.iloc[i-1]

            # Convert to first quadrant since it is symmetric and this cannot occur while passing through an axis
            local1 = [abs(data1['x']), abs(data1['y'])]
            local2 = [abs(data2['x']), abs(data2['y'])]
            global1 = np.abs(data1['global_pos'])
            global2 = np.abs(data2['global_pos'])
            x1, y1 = local1[0] + 2*a*global1[0], local1[1] + 2*b*global1[1]
            x2, y2 = local2[0] + 2*a*global2[0], local2[1] + 2*b*global2[1]

            m = (y2-y1)/(x2-x1)
            y = m*(a-x1) + y1

            if y > b:
                flag = 1
            elif y < b:
                flag = -1
            else:
                flag = 0
            data.at[i, 'symbol'] = symbol + (flag,)
    return data

def batch_run(system, particles, n):
    '''
    Given a system and a set of particles, runs the simulation for each particle for n steps.

    Parameters:
    system (BilliardsSystem): BilliardsSystem object
    particles (list): List of Particle objects
    n (int): Number of steps to run the simulation for each particle

    Returns:
    data (pd.DataFrame): DataFrame containing the data for each particle
    '''

    data = pd.DataFrame()

    for particle in particles:
        system.change_particle(particle)
        d = pd.DataFrame(system.run_simulation(n).iloc[-1]).T
        data = pd.concat([data, d], ignore_index=True)

    return data