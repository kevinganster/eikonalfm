import numpy as np
import pytest
import eikonalfm


@pytest.fixture(scope="module")
def c():
    return np.ones((10, 10))

@pytest.fixture(scope="module")
def x_s():
    return 0, 0

@pytest.fixture(scope="module")
def dx():
    return 0.1, 0.1


@pytest.mark.parametrize("func", [eikonalfm.fast_marching, eikonalfm.factored_fast_marching])
class TestOrder:
    def test_wrong_type(self, func, c, x_s, dx):
        with pytest.raises(TypeError):
            func(c, x_s, dx, None)
        with pytest.raises(TypeError):
            func(c, x_s, dx, 1.2)

    def test_wrong_value(self, func, c, x_s, dx):
        with pytest.raises(ValueError):
            func(c, x_s, dx, 3)

@pytest.mark.parametrize("func", [eikonalfm.fast_marching, eikonalfm.factored_fast_marching])
class TestVelocity:
    def test_no_array(self, func, x_s, dx):
        with pytest.raises(ValueError):
            func(None, x_s, dx, 1)
        with pytest.raises(ValueError):
            func(1.0, x_s, dx, 2)

    def test_dtype(self, func, x_s, dx):
        # string array (int will be cast to float)
        with pytest.raises(ValueError):
            func(np.array([["a", "b"], ["c", "d"]]), x_s, dx, 1)

    def test_dim(self, func, x_s, dx):
        with pytest.raises(ValueError):
            # 1D array, but 2D x_s/dx
            func(np.ones((10, )), x_s, dx, 1)
        with pytest.raises(ValueError):
            # 4D array, but 2D x_s/dx
            func(np.ones((10, 10, 10, 10)), x_s, dx, 2)

@pytest.mark.parametrize("func", [eikonalfm.fast_marching, eikonalfm.factored_fast_marching])
class TestSource:
    def test_no_array(self, func, c, dx):
        with pytest.raises(ValueError):
            func(c, None, dx, 1)
        with pytest.raises(ValueError):
            func(c, 1.0, dx, 2)
    
    def test_dtype(self, func, c, dx):
        with pytest.raises(ValueError):
            func(c, ("a", "b"), dx, 1)
        with pytest.raises(ValueError):
            func(c, np.array([1.5, 2]), dx, 2)
    
    def test_dim(self, func, c, dx):
        with pytest.raises(ValueError):
            func(c, [(0, 0), (1, 1)], dx, 1)
    
    def test_shape(self, func, c, dx):
        # wrong shape
        with pytest.raises(ValueError):
            func(c, (0, ), dx, 1)
        with pytest.raises(ValueError):
            func(c, (0, 0, 0), dx, 1)
        
        # outside shape
        with pytest.raises(ValueError):
            func(c, (-1, 1), dx, 1)
        with pytest.raises(ValueError):
            func(c, (1, -2), dx, 1)
        with pytest.raises(ValueError):
            func(c, (c.shape[0] + 2, 5), dx, 1)
        with pytest.raises(ValueError):
            func(c, (5, c.shape[1] + 2), dx, 1)

@pytest.mark.parametrize("func", [eikonalfm.fast_marching, eikonalfm.factored_fast_marching])
class TestDx:
    def test_no_array(self, func, c, x_s):
        with pytest.raises(ValueError):
            func(c, x_s, None, 1)
        with pytest.raises(ValueError):
            func(c, x_s, 1.3, 2)
    
    def test_dtype(self, func, c, x_s):
        with pytest.raises(ValueError):
            func(c, x_s, ("a", "b"), 1)
    
    def test_dim(self, func, c, x_s):
        with pytest.raises(ValueError):
            func(c, x_s, [(0.1, 0.1), (0.2, 0.2)], 1)
    
    def test_shape(self, func, c, x_s):
        # wrong shape
        with pytest.raises(ValueError):
            func(c, x_s, (0.1, ), 1)
        with pytest.raises(ValueError):
            func(c, x_s, (0.1, 0.1, 0.1), 1)

        # negative/zero
        with pytest.raises(ValueError):
            func(c, (-0.1, 0.1), dx, 1)
        with pytest.raises(ValueError):
            func(c, (0.1, 0), dx, 1)