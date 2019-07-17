import numpy as np
import eikonalfm
import matplotlib.pyplot as plt
from mesh import Mesh2D
print("numpy version:", np.version.version)

np.set_printoptions(precision=2, linewidth=1e10)


mesh = Mesh2D((0, 10), (0, 10), nx1=101)
print(mesh)

# source position
x_s = np.array([0, 0])
# velocity
c = np.ones_like(mesh.X[0])
wall = 1e-8
c[6, :-4] = wall
c[12, :40] = wall
c[12, 41:] = wall
c[22:50, 32] = wall
c[50, 32:80] = wall


tau_fm = eikonalfm.fast_marching(c, x_s, mesh.dx, 2)
tau_ffm = eikonalfm.factored_fast_marching(c, x_s, mesh.dx, 2)
print("finished solving")
print("fm solution:")
print(np.array2string(tau_fm.T, separator="    "))
print("ffm solution:")
print(np.array2string(tau_ffm.T, separator="    "))

tau_fm[tau_fm > 100] = np.nan
tau_ffm[tau_ffm > 100] = np.nan


plt.figure() 
ax = plt.subplot(2, 2, 2)
ax.invert_yaxis()
ax.set_title("Velocity c")
plot = ax.contourf(*mesh.X, c)
plt.colorbar(plot, ax=ax)

  
ax = plt.subplot(2, 2, 3)
ax.invert_yaxis()
ax.set_title("Fast Marching")
# plot = ax.contourf(*mesh.X, tau_fm)
plot = ax.pcolormesh(*mesh.X, tau_ffm, cmap="seismic")
plt.colorbar(plot, ax=ax)

ax = plt.subplot(2, 2, 4)
ax.invert_yaxis()
ax.set_title("Factored Fast Marching")
# plot = ax.contourf(*mesh.X, tau_ffm)
plot = ax.pcolormesh(*mesh.X, tau_ffm, cmap="gist_ncar")
plt.colorbar(plot, ax=ax)
  
plt.show()
