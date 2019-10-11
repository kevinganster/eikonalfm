"""
Eikonal Fast Marching
--------

eikonalfm is a Python (C++-)extension which implements the fast marching method for the eikonal equation
    | grad(tau)(x) |^2 = 1 / c^2(x),
and the factored eikonal equation
    | (tau0 grad(tau1) + tau1 grad(tau0))(x) |^2 = 1 / c^2(x).

See https://github.com/Daarknes/eikonalfm for more information
"""
__version__ = "1.0.0"

from .cfm import fast_marching, factored_fast_marching