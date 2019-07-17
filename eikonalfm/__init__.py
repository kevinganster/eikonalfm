"""
Eikonal Fast Marching
--------

eikonalfm is a Python (C++-)extension which implements the fast marching method for the eikonal equation
    | grad tau(x) |^2 = 1 / c(x)^2.

"""

from .cfm import fast_marching, factored_fast_marching