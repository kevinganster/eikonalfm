import numpy as np
import matplotlib.pyplot as plt
import eikonalfm
from tests.mesh import Mesh2D

np.set_printoptions(precision=2, linewidth=1e10)


mesh = Mesh2D((0, 40), (-40, 0), nx1=101)
print(mesh)

# source position
x_s = np.array([0, 0])
# velocity
c = 2 + np.cos(mesh.X[1])

tau_fm = eikonalfm.fast_marching(c, mesh.pos_to_index(x_s), mesh.dx, 2)
tau_ffm = eikonalfm.factored_fast_marching(c, mesh.pos_to_index(x_s), mesh.dx, 2)


levels = list(range(2, 21, 2))
 
ax = plt.subplot(1, 2, 1)
ax.set_aspect('equal', 'box')
plot = ax.contour(*mesh.X, tau_fm, levels=levels)
ax.clabel(plot, plot.levels, inline_spacing=2, fontsize=12, fmt="%i")
#plt.colorbar(plot, ax=ax)
 
 
ax = plt.subplot(1, 2, 2)
ax.set_aspect('equal', 'box')
plot = ax.contour(*mesh.X, tau_ffm, levels=levels)

    
# start_points = plot.collections[-1].get_paths()[0].vertices[::4]
# ax.streamplot(mesh.X[0].T, mesh.X[1].T, *np.gradient(tau_ffm), start_points=start_points, density=35)
ax.streamplot(mesh.X[0].T, mesh.X[1].T, *np.gradient(tau_ffm), density=1)
#ax.clabel(plot, plot.levels, inline_spacing=12, fontsize=12, fmt="%i")
# ax.plot(*start_points.T, marker="o", markersize=3)
 
 
plt.tight_layout()
# plt.savefig("fig_distance.pdf", bbox_inches="tight", pad_inches=0)
plt.show()