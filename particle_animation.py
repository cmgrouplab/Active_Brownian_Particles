import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

TICK = 10  # Number of steps between each save file
MAX_STEPS = 1000  # Number of movement steps
RADIUS = 0.1  # Radius of particles
LENGTH = 1  # Length of box

INTERVAL = 100  # Time interval between frames in ms
FIGSIZE = 10


def animate_data(data, region, savedir):
    plt.rcParams['animation.html'] = 'html5'

    fig, ax = plt.subplots(figsize=(FIGSIZE, FIGSIZE))

    xmin, xmax, ymin, ymax = region
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal')

    def init():
        return []

    patches = []

    def animate(i):
        if patches:
            for circle, position in zip(patches, data[i]):
                circle.center = position[0], position[1]
        else:
            for position in data[i]:
                patches.append(ax.add_patch(plt.Circle(
                    (position[0], position[1]), RADIUS, edgecolor='b', facecolor='y')))

        ax.set_title(f"Step={i*TICK}")
        return patches

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(data), interval=INTERVAL, blit=True)

    if not os.path.exists(savedir):
        os.mkdir(savedir)

    anim.save(os.path.join(savedir, 'Brownian.gif'),
              dpi=80, writer='imagemagick')


parser = argparse.ArgumentParser()
parser.add_argument('--data-dir', type=str, default='')
parser.add_argument('--save-dir', type=str, default='')
args = parser.parse_args()


files = [f"position{i}.txt" for i in range(0, MAX_STEPS+1, TICK)]
all_data = [np.loadtxt(os.path.join(args.data_dir, f)) for f in files]

animate_data(all_data, [0, LENGTH, 0, LENGTH], args.save_dir)
