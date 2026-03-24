import numpy as np
import eikonalfm
import matplotlib.pyplot as plt


x, dx = np.linspace(-5, 5, 1001, retstep=True)
y, dy = np.linspace(0, 5, 501, retstep=True)
X, Y = np.meshgrid(x, y, indexing="ij")
DX = (dx, dy)
c = 0.5 + 0.1 * Y
x_s = (501, 51)

tau0 = eikonalfm.distance(c.shape, DX, x_s, indexing="ij")
tau1 = eikonalfm.factored_fast_marching(c, x_s, DX, 2)
plt.gca().set_aspect("equal")
mappable = plt.pcolormesh(X, Y, tau0 * tau1)
plt.colorbar(mappable)
# plt.contourf((tau0 * tau1).T)
plt.show()
