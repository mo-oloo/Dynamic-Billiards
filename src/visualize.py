import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import LineCollection
import pandas as pd
import classes as cl

def process_collision_data_incvec(row):
    vec = row['incidence_vector']
    x = row['x']
    y = row['y']

    phi = np.arctan2(vec[1], vec[0]) - np.arctan2(y, x)
    return np.arctan(np.tan(phi))

def process_collision_data(data):
    data = data.dropna(subset=['scatterer_hit'])
    data.loc[:, 'scatterer_hit'] = data['scatterer_hit'].astype(int)
    data.loc[:, 'theta'] = data['theta'].astype(float) % (2*np.pi)

    # Only if incidence_vector isn't already a scalar
    if isinstance(data['incidence_vector'].iloc[0], np.ndarray):
        data.loc[:, 'incidence_vector'] = data.apply(process_collision_data_incvec, axis=1)

    #---------------------------------------
    # Manually combine corner scatterers into a single scatterer
    scat0 = [0]
    scat1 = [1, 2, 3, 4]

    data.loc[data['scatterer_hit'].isin(scat1), 'scatterer_hit'] = 1
    #---------------------------------------

    numScatterers = data['scatterer_hit'].unique()

    return data, numScatterers

def plot_collisions(data, numScatterers, plot_options={}):

    if 'markersize' in plot_options:
        markersize = plot_options['markersize']
    else:
        markersize = 1

    if 's' in plot_options:
        s = plot_options['s']
    else:
        s = None

    fig, axs = plt.subplots(len(numScatterers), 1)
    
    for ax in axs:
        fig_offset = np.pi/8
        ax.set_xlim(0 - fig_offset, 2*np.pi + fig_offset)
        ax.set_ylim(-np.pi/2 - fig_offset, np.pi/2 + fig_offset)
        ax.set_aspect('equal')

    for i in range(len(data)):
        scatterer_hit = data['scatterer_hit'].iloc[i]
        theta = data['theta'].iloc[i]
        incidence_vector = data['incidence_vector'].iloc[i]

        # find the scatterer by using index of numScatterers
        axi = axs[numScatterers.tolist().index(scatterer_hit)]
        axi.plot(theta, incidence_vector, 'ro', markersize=markersize)

        # if i < 2000 and i % 2 == 0:
        #     axi.plot(theta, incidence_vector, 'ro', markersize=markersize)
        # elif i < 2000 and i % 2 == 1:
        #     axi.plot(theta, incidence_vector, 'bo', markersize=markersize)
        # else:
        #     axi.plot(theta, incidence_vector, 'go', markersize=markersize)

    if s is not None:
        fig.suptitle(f'$S_{s}$')

    return fig

def plot_collisions2(data):
    data, numScatterers = process_collision_data(data)
    fig = plot_collisions(data, numScatterers)
    return fig

def plot_trajectories(data):
    pass

def animate_collisions(init_data, data, scatterers, boundary, fps=30, length=10):
    fig, ax = plt.subplots(2, 2, figsize=(6, 6))

    scatterer_circles = [plt.Circle((scatterer.pos_r[0], scatterer.pos_r[1]), scatterer.pos_r[2], color='r', fill=False) for scatterer in scatterers]
    for axi in (ax[0, 0], ax[1, 0]):
        axi.set_xlim(-boundary.width/2, boundary.width/2)
        axi.set_ylim(-boundary.height/2, boundary.height/2)
        # for scatterer_circle in scatterer_circles:
        #     axi.add_patch(scatterer_circle)

    for axi in (ax[0, 1], ax[1, 1]):
        fig_offset = np.pi/8
        axi.set_xlim(0 - fig_offset, 2*np.pi + fig_offset)
        axi.set_ylim(-np.pi/2 - fig_offset, np.pi/2 + fig_offset)
        axi.set_aspect('equal')

    # traj = ax[0, 0].scatter([], [], color='black', s=1)
    # scat_traj_0 = ax[0, 1].scatter([], [], color='black', s=1)
    # scat_traj_1 = ax[1, 1].scatter([], [], color='black', s=1)

    # def init():
    #     traj.set_data([], [])
    #     scat_traj_0.set_data([], [])
    #     scat_traj_1.set_data([], [])
    #     return traj, scat_traj_0, scat_traj_1

    def animate(i):
        ax[0, 0].scatter(init_data['x'].iloc[i], init_data['y'].iloc[i], color='black', s=1)
        ax[1, 0].scatter(data['x'].iloc[i], data['y'].iloc[i], color='black', s=1)
        if data['scatterer_hit'].iloc[i] == 0:
            ax[0, 1].scatter(data['theta'].iloc[i], data['incidence_vector'].iloc[i], color='black', s=1)
        elif data['scatterer_hit'].iloc[i] == 1:
            ax[1, 1].scatter(data['theta'].iloc[i], data['incidence_vector'].iloc[i], color='black', s=1)


    ani = animation.FuncAnimation(fig, animate, frames=int(length*fps), interval=int(1000/fps))
    
    return ani

def animate_trajectories(data, scatterers, boundary, fps=30, length=10):
    fig, ax = plt.subplots()
    ax.set_xlim(-boundary.width/2, boundary.width/2)
    ax.set_ylim(-boundary.height/2, boundary.height/2)
    ax.set_aspect('equal')

    scatterer_circles = [plt.Circle((scatterer.pos_r[0], scatterer.pos_r[1]), scatterer.pos_r[2], color='r', fill=False) for scatterer in scatterers]
    for scatterer_circle in scatterer_circles:
        ax.add_artist(scatterer_circle)

    particle, = ax.plot([], [], 'r', lw=2)

    def init():
        particle.set_data([], [])
        return particle,

    def animate(i):
        particle.set_data(data['x'].iloc[:i], data['y'].iloc[:i])
        if i > 0 and data['scatterer_hit'].iloc[i] is not None and data['scatterer_hit'].iloc[i+1] is not None:
            pass
        else:
            pass
        return particle,

    ani = animation.FuncAnimation(fig, animate, init_func=init, frames=int(length*fps), interval=int(1000/fps), blit=True)
    
    return ani