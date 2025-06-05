import numpy as np

import eikonalfm
np.set_printoptions(linewidth=int(1e6))


def test_fm():
    c = 0.5 * np.ones((101, 101))
    x_s = 5, 5
    dx = 0.1, 0.2

    tau = eikonalfm.fast_marching(c, x_s, dx, 2)
    assert isinstance(tau, np.ndarray) and tau.dtype == float and tau.shape == c.shape
    # can't really test tau values since the results of non-factored FM are not good enough

    tau, sequence, orders = eikonalfm.fast_marching(c, x_s, dx, 2, output_sensitivities=True)
    assert isinstance(tau, np.ndarray) and tau.dtype == float and tau.shape == c.shape
    assert isinstance(sequence, np.ndarray) and sequence.dtype == int and sequence.shape == c.shape
    assert isinstance(orders, np.ndarray) and orders.dtype == np.int8 and orders.shape == (c.ndim, *c.shape)

def test_ffm():
    c = 0.5 * np.ones((101, 101))
    x_s = 5, 5
    dx = 0.1, 0.2

    tau1 = eikonalfm.factored_fast_marching(c, x_s, dx, 2)
    assert isinstance(tau1, np.ndarray) and tau1.dtype == float and tau1.shape == c.shape

    tau1, sequence, orders = eikonalfm.fast_marching(c, x_s, dx, 2, output_sensitivities=True)
    assert isinstance(tau1, np.ndarray) and tau1.dtype == float and tau1.shape == c.shape
    assert isinstance(sequence, np.ndarray) and sequence.dtype == int and sequence.shape == c.shape
    assert isinstance(orders, np.ndarray) and orders.dtype == np.int8 and orders.shape == (c.ndim, *c.shape)