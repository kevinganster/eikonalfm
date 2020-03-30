"""
Eikonal Fast Marching
--------

eikonalfm is a Python (C++-)extension which implements the fast marching method for the eikonal equation
    | grad(tau)(x) |^2 = 1 / c^2(x),
and the factored eikonal equation
    | (tau0 grad(tau1) + tau1 grad(tau0))(x) |^2 = 1 / c^2(x).

See https://github.com/Daarknes/eikonalfm for more information
"""
__version__ = "0.9.3"

import numpy as np
from .cfm import fast_marching, factored_fast_marching


def distance(shape, dx, x_s, indexing="ij"):
    """
    Calculates a distance field according to the euclidian metric.

    Parameters
    ----------
        shape : sequence of ints
            Shape of the generated array.
        x_s : sequence of doubles
            Source position as position-vector, e.g. ``(0.5, 1)``.
            Must have the same length as the number of dimensions of shape.
        dx : sequence of doubles
            Grid spacing for each dimension, dx > 0.
            Must have the same length as the number of dimensions of shape.

    Returns
    ----------
        tau0 : ndarray
            euclidian distance field |x - x_s|.
    """
    x = []
    for shape_i, dx_i in zip(shape, dx):
        x.append(np.linspace(0, shape_i*dx_i, shape_i, endpoint=False))

    mesh = np.array(np.meshgrid(*x, indexing=indexing))
    return np.linalg.norm((mesh.T - x_s).T, ord=2, axis=0)
