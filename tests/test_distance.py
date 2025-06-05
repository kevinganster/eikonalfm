import itertools
import numpy as np
import eikonalfm
# np.set_printoptions(linewidth=int(1e6))


def test_dist():
    shape = 101, 101
    dx = 0.1, 0.2
    x_s = np.array((50, 50))
    d = eikonalfm.distance(shape, dx, x_s, indexing="ij")

    assert isinstance(d, np.ndarray)
    assert np.allclose(eikonalfm.distance(shape, dx, x_s, indexing="xy"), d.T)

    # check values
    for i in itertools.product(*(range(n) for n in shape)):
        assert np.isclose(np.linalg.norm((i - x_s) * dx), d[i])