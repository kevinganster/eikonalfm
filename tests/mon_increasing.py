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


# procedures.test_single(((0, 8), (0, 4)), 321, c_func, tau_exact, x_s=(4, 0), order=2)
# procedures.test_isochrones(((0, 8), (0, 4)), 321, c_func, tau_exact, x_s=(3, 0), x_r=(5, 0), order=2)
# procedures.test_convergence(((0, 8), (0, 4)), (11, 21, 41, 81, 161, 321, 641, 1281, 2561), c_func, tau_exact, x_s=(4, 0), order=2)
procedures.test_convergence(((0, 8), (0, 4)), (321, 641, 1281, 2561, 5121, 10241), c_func, tau_exact, x_s=(4, 0), order=2)
