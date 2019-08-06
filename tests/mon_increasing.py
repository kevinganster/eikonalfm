import numpy as np
from tests import procedures

np.set_printoptions(precision=2, linewidth=1e10)


# non-convergent solution (decreasing velocity):
# a = 0.2
# b = 2
# mesh = Mesh2D((0, 18), (0, 25), nx1=181)

# Haber's constant gradient of squared slowness
a = -0.4
b = 2


def tau_exact(x, y):
    S = b**2 + a * np.abs(x[1] - y[1])
    norm_sq = (x[0] - y[0])**2 + (x[1] - y[1])**2
    sigma = (2 * norm_sq) / (S + np.sqrt(S**2 - a**2 * norm_sq))
    return S * np.sqrt(sigma) -  a**2 / 6 * sigma**(3 / 2)

def c_func(x):
    return 1 / np.sqrt(b**2 + 2*a * x[1])


area = ((0, 8), (0, 4))
# procedures.test_single(area, 321, c_func, tau_exact, x_s=(4, 0), order=2)
# procedures.test_isochrones(area, 321, c_func, tau_exact, x_s=(3, 0), x_r=(5, 0), order=2)
# procedures.test_convergence(area, (11, 21, 41, 81, 161, 321, 641, 1281, 2561), c_func, tau_exact, x_s=(4, 0), order=2)
# procedures.test_convergence(area, (321, 641, 1281, 2561, 5121, 10241), c_func, tau_exact, x_s=(4, 0), order=2)


from tests.mesh import Mesh2D
import eikonalfm
#np.random.seed(42)

order = 2
factor = 2
mesh = mesh = Mesh2D(*area, nx1=321)
x_s = x_s = mesh.pos_to_index(np.array((4, 0)))

c = c_func(mesh.X)
c_ = c + factor * np.random.random(size=c.shape)

#tau = tau_exact(mesh.X, mesh.index_to_pos(x_s))
 
print("c-inf-error:", np.linalg.norm((c - c_).flatten(), ord=np.inf))

tau = eikonalfm.fast_marching(c, x_s, mesh.dx, order)
tau_fm = eikonalfm.fast_marching(c_, x_s, mesh.dx, order)
print("tau_fm-inf-error:", np.linalg.norm((tau - tau_fm).flatten(), ord=np.inf))

tau = eikonalfm.factored_fast_marching(c, x_s, mesh.dx, order)
tau_ffm = eikonalfm.factored_fast_marching(c_, x_s, mesh.dx, order)
print("tau_ffm-inf-error:", np.linalg.norm((tau - tau_ffm).flatten(), ord=np.inf))