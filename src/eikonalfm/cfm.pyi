from typing import TYPE_CHECKING, Any, Literal, Sequence, Union
from numpy.typing import ArrayLike, NDArray

if TYPE_CHECKING or sys.version_info >= (3, 9):
    import numpy as np
    IntArray = NDArray[np.int_]
    FloatArray = NDArray[np.floating[Any]]
else:
    IntArray = NDArray
    FloatArray = NDArray


def fast_marching(c: ArrayLike, x_s: Union[Sequence[int], IntArray], dx: Union[Sequence[float], FloatArray], order: Literal[1, 2]) -> FloatArray:
    """
    Calculates amplitude factor of the transport equation.

    Parameters
    ----------
        c : ndarray
            backgound velocity, c(x) > 0.
        x_s : sequence of ints
            Source position as 2D index vector, e.g. ``(0, 0)``.
            Must have the same length as the number of dimensions of c.
        dx : sequence of doubles
            Grid spacing for each dimension, dx > 0.
            Must have the same length as the number of dimensions of c.
        order : {1, 2}
            Order of the finite difference operators.

    Returns
    ----------
        tau : ndarray
            numerical solution tau for the eikonal equation."""
    ...

def factored_fast_marching(c: ArrayLike, x_s: Union[Sequence[int], IntArray], dx: Union[Sequence[float], FloatArray], order: Literal[1, 2]) -> FloatArray:
    """
    Calculates amplitude factor of the transport equation.

    Parameters
    ----------
        c : ndarray
            backgound velocity, c(x) > 0.
        x_s : sequence of ints
            Source position as 2D index vector, e.g. ``(0, 0)``.
            Must have the same length as the number of dimensions of c.
        dx : sequence of doubles
            Grid spacing for each dimension, dx > 0.
            Must have the same length as the number of dimensions of c.
        order : {1, 2}
            Order of the finite difference operators.

    Returns
    ----------
        tau1 : ndarray
            numerical solution tau for the factored eikonal equation."""
    ...
