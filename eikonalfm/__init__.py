"""
Eikonal Fast Marching
--------

eikonalfm is a Python (C++-)extension which implements the fast marching method for the eikonal equation
    | grad(tau)(x) |^2 = 1 / c^2(x),
and the factored eikonal equation
    | (tau0 grad(tau1) + tau1 grad(tau0))(x) |^2 = 1 / c^2(x).

See https://github.com/Daarknes/eikonalfm for more information
"""
__version__ = "0.9.2"

import numpy as np
from .cfm import fast_marching, factored_fast_marching


def distance(shape, dx, x_s, indexing="ij"):
    x = []
    for shape_i, dx_i in zip(shape, dx):
        x.append(np.linspace(0, shape_i*dx_i, shape_i, endpoint=False))

    mesh = np.array(np.meshgrid(*x, indexing=indexing))
    return np.linalg.norm((mesh.T - x_s).T, ord=2, axis=0)
