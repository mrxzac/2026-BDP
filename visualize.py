import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from config_0 import beam_0,u

def plot_matboard(
    board,ax,
    color="skyblue",
    edgecolor="navy",
    alpha=0.4,
    show_labels=False,
):
    vertices = board.vertices
    axis = {"x": 0, "y": 1, "z": 2}[board.dirc]
    offset = np.zeros(3)
    offset[axis] = board.t.to(u.mm).magnitude / 2
    front = vertices + offset
    back = vertices - offset
    n = len(vertices)
    faces = [front, back[::-1]]
    for i in range(n):
        j = (i + 1) % n
        faces.append([front[i],front[j],back[j],back[i],])
    solid = Poly3DCollection(faces,facecolors=color,edgecolors=edgecolor,linewidths=0.7,alpha=alpha)
    ax.add_collection3d(solid)
    return np.vstack([front, back])


def plot_beam(beam,alpha=0.4,padding=5):
    fig = plt.figure(figsize=(15, 7))
    ax = fig.add_subplot(111, projection="3d")
    colors = {"x": "orange","y": "skyblue","z": "lightgreen"}
    all_vertices = []
    for board in beam.matboards:
        vertices = plot_matboard(board,ax,color=colors[board.dirc],alpha=alpha,show_labels=show_labels,)
        all_vertices.append(vertices)

    vertices = np.vstack(all_vertices)
    mins = vertices.min(axis=0) - padding
    maxs = vertices.max(axis=0) + padding
    ax.view_init(elev=25, azim=-60, vertical_axis="y")
    ax.set_xlim(mins[0], maxs[0])
    ax.set_ylim(mins[1], maxs[1])
    ax.set_zlim(mins[2], maxs[2])
    ax.set_aspect("equal", adjustable="box")
    # ax.set_xlabel("X (mm)")
    # ax.set_ylabel("Y (mm)")
    # ax.set_zlabel("Z (mm)")
    # plt.tight_layout()
    plt.show()
    return ax

plot_beam(beam_0, alpha=0.15)