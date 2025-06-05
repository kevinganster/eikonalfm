from typing import TYPE_CHECKING, Any, Literal, Sequence, Tuple, Union, overload
from numpy.typing import ArrayLike, NDArray

if TYPE_CHECKING or sys.version_info >= (3, 9):
    import numpy as np
    IntArray = NDArray[np.int_]
    FloatArray = NDArray[np.floating[Any]]
else:
    IntArray = NDArray
    FloatArray = NDArray


@overload
def fast_marching(
    c: ArrayLike,
    x_s: Union[Sequence[int], IntArray],
    dx: Union[Sequence[float], FloatArray],
    order: Literal[1, 2],
    output_sensitivities: Literal[False] = False
) -> FloatArray:
    r"""
    Calculates the fast marching solution to the eikonal equation.

    Parameters
    ----------
        c : ndarray
            backgound velocity, `c(x) > 0`.
        x_s : sequence of ints
            Source position as index vector, e.g. `(0, 0)`.
            Must have the same length as the number of dimensions of c.
        dx : sequence of doubles
            Grid spacing for each dimension, `dx > 0`.
            Must have the same length as the number of dimensions of c.
        order : {1, 2}
            Order of the finite difference operators.
        output_sensitivities : boolean, optional
            Additionally returns sensitivity data. Default is `False`.

    Returns
    ----------
        tau : ndarray
            numerical solution tau for the eikonal equation.
        sequence : ndarray, optional
            Only returned when `output_sensitivities` is set to `True`. The sequence in which each gridpoint was set to known.
        orders : ndarray, optional
            Only returned when `output_sensitivities` is set to `True`. The finite difference orders for each dimension used at each gridpoint.
    """
    ...
@overload
def fast_marching(
    c: ArrayLike,
    x_s: Union[Sequence[int], IntArray],
    dx: Union[Sequence[float], FloatArray],
    order: Literal[1, 2],
    output_sensitivities: Literal[True]
) -> Tuple[FloatArray, IntArray, IntArray]:
    """
    Calculates the fast marching solution to the eikonal equation.

    Parameters
    ----------
        c : ndarray
            backgound velocity, `c(x) > 0`.
        x_s : sequence of ints
            Source position as index vector, e.g. `(0, 0)`.
            Must have the same length as the number of dimensions of c.
        dx : sequence of doubles
            Grid spacing for each dimension, `dx > 0`.
            Must have the same length as the number of dimensions of c.
        order : {1, 2}
            Order of the finite difference operators.
        output_sensitivities : boolean, optional
            Additionally returns sensitivity data. Default is `False`.

    Returns
    ----------
        tau : ndarray
            numerical solution tau for the eikonal equation.
        sequence : ndarray
            The sequence in which each gridpoint was set to known.
        orders : ndarray
            The finite difference orders for each dimension used at each gridpoint.
    """
    ...

@overload
def factored_fast_marching(
    c: ArrayLike,
    x_s: Union[Sequence[int], IntArray],
    dx: Union[Sequence[float], FloatArray],
    order: Literal[1, 2],
    output_sensitivities: Literal[False] = False
) -> FloatArray:
    """
    Calculates the fast marching solution to the factored eikonal equation.

    Parameters
    ----------
        c : ndarray
            backgound velocity, `c(x) > 0`.
        x_s : sequence of ints
            Source position as index vector, e.g. `(0, 0)`.
            Must have the same length as the number of dimensions of c.
        dx : sequence of doubles
            Grid spacing for each dimension, `dx > 0`.
            Must have the same length as the number of dimensions of c.
        order : {1, 2}
            Order of the finite difference operators.
        output_sensitivities : boolean, optional
            Additionally returns sensitivity data. Default is `False`.

    Returns
    ----------
        tau1 : ndarray
            numerical solution tau1 for the factored eikonal equation.
        sequence : ndarray, optional
            Only returned when `output_sensitivities` is set to `True`. The sequence in which each gridpoint was set to known.
        orders : ndarray, optional
            Only returned when `output_sensitivities` is set to `True`. The finite difference orders for each dimension used at each gridpoint.
    """
    ...

@overload
def factored_fast_marching(
    c: ArrayLike,
    x_s: Union[Sequence[int], IntArray],
    dx: Union[Sequence[float], FloatArray],
    order: Literal[1, 2], output_sensitivities: Literal[True]) -> Tuple[FloatArray, IntArray, IntArray]:
    """
    Calculates the fast marching solution to the eikonal equation.

    Parameters
    ----------
        c : ndarray
            backgound velocity, `c(x) > 0`.
        x_s : sequence of ints
            Source position as index vector, e.g. `(0, 0)`.
            Must have the same length as the number of dimensions of c.
        dx : sequence of doubles
            Grid spacing for each dimension, `dx > 0`.
            Must have the same length as the number of dimensions of c.
        order : {1, 2}
            Order of the finite difference operators.
        output_sensitivities : boolean, optional
            Additionally returns sensitivity data. Default is `False`.

    Returns
    ----------
        tau1 : ndarray
            numerical solution tau1 for the factored eikonal equation.
        sequence : ndarray
            The sequence in which each gridpoint was set to known.
        orders : ndarray
            The finite difference orders for each dimension used at each gridpoint.
    """
    ...
