import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import LineCollection
import pandas as pd
import classes as cl

def plot_collisions(data):
    data = data.dropna(subset=['scatterer_hit'])
    data.loc[:, 'scatterer_hit'] = data['scatterer_hit'].astype(int)
    data.loc[:, 'theta'] = data['theta'].astype(float)
    data.loc[:, 'incidence_vector'] = data['incidence_vector'].apply(lambda x: np.arctan2(x[1], x[0]))

    #---------------------------------------
    # Manually combine corner scatterers into a single scatterer
    scat0 = [0]
    scat1 = [1, 2, 3, 4]

    data.loc[data['scatterer_hit'].isin(scat1), 'scatterer_hit'] = 1
    #---------------------------------------

    numScatterers = data['scatterer_hit'].unique()
    fig, axs = plt.subplots(len(numScatterers), 1)
    
    for ax in axs:
        fig_offset = np.pi/4
        ax.set_xlim(-np.pi - (fig_offset), np.pi + (fig_offset))
        ax.set_ylim(-np.pi/2 - (fig_offset/2), np.pi/2 + (fig_offset/2))
        ax.set_aspect('equal')

    for i in range(len(data)):
        scatterer_hit = data['scatterer_hit'].iloc[i]
        theta = data['theta'].iloc[i]
        incidence_vector = data['incidence_vector'].iloc[i]

        # find the scatterer by using index of numScatterers
        ax = axs[numScatterers.tolist().index(scatterer_hit)]
        ax.plot(incidence_vector, theta, 'ro', markersize=1)

    return fig

def plot_trajectories(data):
    pass

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