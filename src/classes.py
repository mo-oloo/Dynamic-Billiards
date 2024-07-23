import numpy as np
import pandas as pd
import functions as fn
import visualize as vis

class Particle:
    def __init__(self, x=None, y=None, vx=None, vy=None, arr=None):
        if arr is not None:
            self.pos = np.array(arr[:2])
            self.vel = np.array(arr[2:])
        elif x is None and y is None and vx is None and vy is None and arr is None:
            # Creates an 'empty' particle
            self.pos = self.vel = np.empty(2)
        else:
            if x is None or y is None or vx is None or vy is None:
                raise ValueError("Please provide x, y, vx, and vy values.")
            self.pos = np.array([x, y], dtype=np.float64)
            v = np.array([vx, vy], dtype=np.float64)
            self.vel = v / np.linalg.norm(v)

    def update_position(self, t):
        self.pos += t * self.vel
        # self.pos = np.add(self.pos, t * self.vel, casting='unsafe')

    def update_velocity(self, scatterer_pos):
        xp, yp = self.pos
        xc, yc = scatterer_pos
        incidence_vector = self.vel / np.linalg.norm(self.vel) # Normalize velocity vector
        theta = np.arctan2(yp - yc, xp - xc)
        normal_vector = np.array([np.cos(theta), np.sin(theta)])
        reflection_vector = incidence_vector - 2 * np.dot(incidence_vector, normal_vector) * normal_vector
        reflection_vector = reflection_vector / np.linalg.norm(reflection_vector) # Normalize reflection vector
        self.vel = reflection_vector
        phi = np.arctan2(reflection_vector[1], reflection_vector[0]) - theta

        return theta, np.arctan(np.tan(phi))

class Scatterer:
    def __init__(self, x, y, r):
        self.pos = np.array([x, y])
        self.radius = r
        self.pos_r = np.array([x, y, r])

class Boundary:
    def __init__(self, w, h, center=np.array([0, 0])):
        self.origin = center
        self.width = w
        self.height = h
        self.dimensions = np.array([w, h])

class BilliardsSystem:

    def __init__(self, particle, scatterers, boundary, a, b):
        self.particle = particle
        self.scatterers = scatterers
        self.scatpos = np.array([a, b])
        self.boundary = boundary
        self.time = 0
        self.num_collisions = 0
        self.global_pos = np.array([0, 0])
        self.last_collision = np.array([0, 0])
        self.start_data = self.init_start_data()
        self.data = [self.start_data]

    def change_scatterers(self, r1, r2):
        a, b = self.scatpos
        scatterer0 = Scatterer(0, 0, r1)
        scatterer1 = Scatterer(a, b, r2)
        scatterer2 = Scatterer(a, -b, r2)
        scatterer3 = Scatterer(-a, b, r2)
        scatterer4 = Scatterer(-a, -b, r2)
        self.scatterers = [scatterer0, scatterer1, scatterer2, scatterer3, scatterer4]

    def init_start_data(self, scatterer_hit=None, theta=None, reflection_vector=None):
        # IMPORTANT: This function actually uses the current location of the particle, only use when initializing the system
        start_data = {
            'time': self.time,
            'n': self.num_collisions,
            'x': self.particle.pos[0],
            'y': self.particle.pos[1],
            'vx': self.particle.vel[0],
            'vy': self.particle.vel[1],
            'scatterer_hit': scatterer_hit.astype(int) if scatterer_hit is not None else None,
            'theta': theta,
            'reflection_vector': reflection_vector,
            'global_pos': np.array([0, 0]), # Assumes it starts at the origin
            'last_collision': np.array([0, 0]),
            'symbol': np.array([0, 0])
        }
        return start_data

    def reset_data(self):
        self.time = 0
        self.num_collisions = 0
        self.particle.pos = np.array([self.start_data['x'], self.start_data['y']], dtype=np.float64)
        self.particle.vel = np.array([self.start_data['vx'], self.start_data['vy']], dtype=np.float64)
        self.global_pos = np.array([0, 0])
        self.last_collision = np.array([0, 0])
        self.data = [self.start_data]

    def run_simulation(self, n):
        self.reset_data()
        i = 0
        while i < n:
            i += self.update()
        return pd.DataFrame(self.data)
    
    def change_particle(self, particle, scatterer_hit=None, theta=None, reflection_vector=None):
        self.particle = particle
        self.start_data = self.init_start_data(scatterer_hit, theta, reflection_vector)
        self.reset_data()

    def update(self):
        time, scatterer_index = self.find_next_collision()
        self.particle.update_position(time)
        self.time += time

        new_pos, change_vector = fn.check_particle_bounds(self.particle.pos, self.boundary.dimensions)
        if new_pos is not None:
            new_data = self.data_entry()
            self.data.append(new_data)
            self.particle.pos = new_pos

        self.global_pos = self.global_pos + change_vector if change_vector is not None else self.global_pos

        if scatterer_index is None:
            theta, incidence_vector, symbol = None, None, None
            collide_true = 0
        else:
            scatterer = self.scatterers[scatterer_index]
            theta, incidence_vector = self.particle.update_velocity(scatterer.pos)
            collide_true = 1
            self.num_collisions += 1

            # position of scatterer hit
            collision_pos_global = scatterer.pos + 2 * self.global_pos * self.scatpos
            symbol = collision_pos_global - self.last_collision
            self.last_collision = collision_pos_global
        
        new_data = self.data_entry(scatterer_index, theta, incidence_vector, symbol)
        self.data.append(new_data)
        return collide_true

    def find_next_collision(self):
        t = []
        time_collision_boundary = fn.find_collisions_boundary(self.particle.pos, self.particle.vel, self.boundary.dimensions)
        if time_collision_boundary is not None:
            t.append((time_collision_boundary, None))

        for scatterer in self.scatterers:
            ti = fn.find_collision_scatterer(self.particle.pos, self.particle.vel, scatterer.pos_r)
            if ti is not None:
                t.append((ti, self.scatterers.index(scatterer)))

        if t:
            t_min, i_min = min(t, key=lambda x: x[0])
            return t_min, i_min
        else: # Shouldn't happen, since particle will always collide with boundary
            SystemError("No collision found")
        
    def data_entry(self, scat_index=None, theta=None, reflection_vector=None, symbol=None):
        new_data = {
                'time': self.time,
                'n': self.num_collisions,
                'x': self.particle.pos[0],
                'y': self.particle.pos[1],
                'vx': self.particle.vel[0],
                'vy': self.particle.vel[1],
                'scatterer_hit': scat_index,
                'theta': theta,
                'reflection_vector': reflection_vector,
                'global_pos': self.global_pos,
                'last_collision': self.last_collision,
                'symbol': symbol
        }
        return new_data
        
    # Graph angle of incidence vs theta by scatterer
    def plot_collisions(self):
        vis.plot_collisions(self.data)

    # Visualize trajectories
    def plot_trajectories(self):
        vis.plot_trajectories(self.data)

    # Animate trajectories
    def animate_trajectories(self, fps=30, length=None):
        if length is None:
            length = int(np.ceil(self.data['time'].iloc[-1]))
        return vis.animate_trajectories(self.data, self.scatterers, self.boundary, fps, length)

    
## ADDITIONAL FUNCTIONS
    
def random_particle(scatterers, boundary):
    '''
    Generates a particle with a random position and initial velocity within the boundary and
    outside of the scatterers.
    '''
    while True:
        x = np.random.uniform(-boundary.width/2, boundary.width/2)
        y = np.random.uniform(-boundary.height/2, boundary.height/2)
        if all(np.linalg.norm(np.array([x, y]) - scatterer.pos) > scatterer.radius for scatterer in scatterers):
            break
    vx = np.random.uniform(-1, 1)
    vy = np.random.uniform(-1, 1)
    return Particle(x, y, vx, vy)

