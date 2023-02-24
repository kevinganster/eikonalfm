"""
Eikonal Fast Marching
--------

eikonalfm is a Python (C++-)extension which implements the fast marching method for the eikonal equation
    | grad(tau)(x) |^2 = 1 / c^2(x),
and the factored eikonal equation
    | (tau0 grad(tau1) + tau1 grad(tau0))(x) |^2 = 1 / c^2(x).

See https://github.com/kevinganster/eikonalfm for more information
"""
__title__ = "eikonalfm"
__version__ = "0.9.7"
__author__ = "Kevin Ganster"
__email__ = "kevinganster@gmail.com"


import numpy as np
from .cfm import fast_marching, factored_fast_marching


def distance(shape, dx, x_s, indexing="xy"):
    """
    Calculates a distance field according to the euclidian metric.

    Parameters
    ----------
        shape : sequence of ints
            Shape of the generated array.
        x_s : sequence of ints
            Source position as index-vector relative to the shape, e.g. ``(9, 3)``.
            Must have the same length as the number of dimensions of shape.
        dx : sequence of doubles
            Grid spacing for each dimension, dx > 0.
            Must have the same length as the number of dimensions of shape.
        indexing : {'xy', 'ij'}, optional
            Cartesian ('xy', default) or matrix ('ij') indexing of output.
            See numpy's meshgrid function for more details.

    Returns
    ----------
        tau0 : ndarray
            euclidian distance field |x - x_s|.
    """
    assert(len(shape) == len(dx) == len(x_s))
    x = []
    for shape_i, dx_i in zip(shape, dx):
        x.append(np.linspace(0, shape_i*dx_i, shape_i, endpoint=False))
    mesh = np.array(np.meshgrid(*x, indexing=indexing))
    
    # get real coordinates of x_s
    x_s = np.array(x_s) * np.array(dx)
    return np.linalg.norm((mesh.T - x_s).T, ord=2, axis=0)
