import numpy as np
import matplotlib.pyplot as plt
import eikonalfm
from tests.mesh import Mesh2D


def tau_exact(x, y):
    return ((x[0] - y[0])**2 + (x[1] - y[1])**2)**0.5

mesh = Mesh2D((-4, 4), (-4, 4), nx1=41)
print(mesh)

# source position
x_s = np.array([0, 0])
# velocity
c = np.ones_like(mesh.X[0])

tau = tau_exact(mesh.X, x_s)
tau_fm = eikonalfm.fast_marching(c, mesh.pos_to_index(x_s), mesh.dx, 2)
tau_ffm = eikonalfm.factored_fast_marching(c, mesh.pos_to_index(x_s), mesh.dx, 2)


levels = list(range(1, 7))
ax = plt.gca()
ax.set_aspect('equal', 'box')

#plot = ax.contourf(*mesh.X, tau, alpha=0.2, levels=[0]+levels)
#plot = ax.pcolormesh(*mesh.X, tau, cmap="gray", alpha=0.3, shading="gouraud")
#plt.colorbar(plot, ax=ax)

plot = ax.contour(*mesh.X, tau, levels=levels)
ax.clabel(plot, plot.levels, inline_spacing=14, fontsize=14, fmt="%i")

ax.plot(*x_s, color="red", marker="o", markersize=8)

# remove axis labels
ax.set_xticks([])
ax.set_yticks([])

plt.tight_layout()

# plt.savefig("fig_distance.pdf", bbox_inches="tight", pad_inches=0)
plt.show()