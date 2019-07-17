import numpy as np
from tests import procedures

np.set_printoptions(precision=2, linewidth=1e10)


def tau_exact(x, y):
    return np.full_like(x[0], np.nan, dtype=np.double)

def c_func(x):
    return 1 + x[1] / 8 * (np.sin(8*np.pi*x[1]) + np.cos(6*np.pi*x[0]))


area = ((0, 1), (0, 1))
procedures.test_single(area, 1001, c_func, tau_exact, x_s=(0, 0), order=2)
procedures.test_isochrones(area, 1001, c_func, tau_exact, x_s=(0, 0), x_r=(1, 0), order=2)
# procedures.test_convergence(area, (11, 21, 41, 81, 161, 321, 641, 1281, 2561), c_func, tau_exact, order=2)
# procedures.test_convergence(area, (321, 641, 1281, 2561, 5121, 10241), c_func, tau_exact, order=2)
