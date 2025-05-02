import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
from numpy import genfromtxt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

def visualizeCurrentCircularGeometry(radius: float,
                                     deviceData,
                                     ): 
    
    
    
    hopping_counts = deviceData["event_counts"]
    acceptor_coords = deviceData["acceptor_coordinates"]
    donor_coords = deviceData["donor_coordinates"]
    total_time = deviceData["device_time"]
    
    net_hops = (hopping_counts - hopping_counts.T) / total_time

    lower_half_vals = np.abs(np.tril(net_hops).flatten())
    min_current = lower_half_vals.min()
    max_current = lower_half_vals.max()

    current_range = np.linspace(0.0, max_current, 256)
    colors = np.array([[0, 0, 0, np.sqrt(val/max_current)] for val in current_range])
    custom_map = ListedColormap(colors)
    norm = mcolors.Normalize(vmin=min_current, vmax=max_current)
    sm = ScalarMappable(norm=norm, cmap=custom_map)
    sm.set_array([])

    padding = 0.1*radius
    x_lim = (-radius - padding, radius + padding)
    y_lim = (-radius - padding, radius + padding)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(2)

    ax.scatter(acceptor_coords[:, 0], acceptor_coords[:, 1], c="k", zorder=3)    
    ax.scatter(donor_coords[:, 0], donor_coords[:, 1], c="dodgerblue", zorder=3) 

    device_boundary = patches.Circle((0,0), radius=radius+0.05, fill=False, edgecolor="gray", lw=5)
    ax.add_patch(device_boundary)

    n_acceptors = acceptor_coords.shape[0]
    for i in range(n_acceptors):
        for j in range(i+1, n_acceptors):
            current_value = np.abs(net_hops[i, j])
            if current_value > 0:
                ax.plot([acceptor_coords[i, 0], acceptor_coords[j, 0]],
                        [acceptor_coords[i, 1], acceptor_coords[j, 1]],
                        color="black", lw=2,
                        alpha=np.sqrt(current_value / max_current))

    plt.show()

if __name__ == "__main__":
    radius = 150.0
    R = np.sqrt(np.pi*radius*radius / 200)
    radius = radius / R
    data = np.load("currentData/device_1.npz")
    visualizeCurrentCircularGeometry(radius, data)