import numpy as np
import pandas as pd
import functions as fn
import visualize as vis

class Particle:
    def __init__(self, x, y, vx, vy):
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
        normal_vector = np.array([-np.cos(theta), -np.sin(theta)])
        reflection_vector = incidence_vector - 2 * np.dot(self.vel, normal_vector) * normal_vector
        reflection_vector = reflection_vector / np.linalg.norm(reflection_vector) # Normalize reflection vector
        self.vel = reflection_vector

        return theta, incidence_vector

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

    def __init__(self, particle, scatterers, boundary):
        self.particle = particle
        self.scatterers = scatterers
        self.boundary = boundary
        self.time = 0
        self.num_collisions = 0
        self.start_data = {
            'time': [self.time],
            'n': [self.num_collisions],
            'x': [self.particle.pos[0]],
            'y': [self.particle.pos[1]],
            'vx': [self.particle.vel[0]],
            'vy': [self.particle.vel[1]],
            'scatterer_hit': [None],
            'theta': [None],
            'incidence_vector': [None]
        }
        self.data = pd.DataFrame(data=self.start_data)

    def reset_data(self):
        self.time = 0
        self.num_collisions = 0
        self.particle.pos = np.array([self.start_data['x'][0], self.start_data['y'][0]], dtype=np.float64)
        self.particle.vel = np.array([self.start_data['vx'][0], self.start_data['vy'][0]], dtype=np.float64)
        self.data = pd.DataFrame(data=self.start_data)


    def run_simulation(self, n):
        self.reset_data()
        i = 0
        while i < n:
            i += self.update()
        return self.data

    def update(self):
        time, scatterer_index = self.find_next_collision()
        self.particle.update_position(time)
        self.time += time

        new_pos = fn.check_particle_bounds(self.particle.pos, self.boundary.dimensions)
        if new_pos is not None:
            new_data = self.data_entry()
            self.data = pd.concat([self.data, new_data], ignore_index=True)
            self.particle.pos = new_pos

        if scatterer_index is None:
            theta, incidence_vector = None, None
            collide_true = 0
        else:
            scatterer = self.scatterers[scatterer_index]
            theta, incidence_vector = self.particle.update_velocity(scatterer.pos)
            collide_true = 1
            self.num_collisions += 1
        
        new_data = self.data_entry(scatterer_index, theta, incidence_vector)
        self.data = pd.concat([self.data, new_data], ignore_index=True)
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
            print("THERE WAS AN ERROR")
            return None, None
        
    def data_entry(self, scat_index=None, theta=None, incidence_vector=None):
        new_data = {
                'time': [self.time],
                'n': [self.num_collisions],
                'x': [self.particle.pos[0]],
                'y': [self.particle.pos[1]],
                'vx': [self.particle.vel[0]],
                'vy': [self.particle.vel[1]],
                'scatterer_hit': [scat_index],
                'theta': [theta],
                'incidence_vector': [incidence_vector]
        }
        new_df = pd.DataFrame(data=new_data)
        return new_df
        
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

