import numpy as np
from tests import procedures

np.set_printoptions(precision=2, linewidth=1e10)


# Haber's constant gradient of velocity
a = 1.0
b = 0.5

def tau_exact(x, y):
    norm_sq = ((x[0] - y[0])**2 + (x[1] - y[1])**2)
    return 1 / a * np.arccosh(1 + a**2 / (2*b) * norm_sq / (b + a * np.abs(x[1] - y[1])))

def c_func(x):
    return b + a * x[1]


# procedures.test_single(((0, 8), (0, 4)), 321, c_func, tau_exact, x_s=(4, 0), order=2)
# procedures.test_isochrones(((0, 8), (0, 4)), 321, c_func, tau_exact, x_s=(3, 0), x_r=(5, 0), order=2)
# procedures.test_convergence(((0, 8), (0, 4)), (11, 21, 41, 81, 161, 321, 641, 1281, 2561), c_func, tau_exact, x_s=(4, 0), order=2)
procedures.test_convergence(((0, 8), (0, 4)), (321, 641, 1281, 2561, 5121, 10241), c_func, tau_exact, x_s=(4, 0), order=2)
