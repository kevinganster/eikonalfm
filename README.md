# Eikonal Fast Marching

eikonalfm is a Python (C++) extension which implements the Fast Marching method for the eikonal equation  
<p align="center"><img src="https://latex.codecogs.com/svg.latex?\lvert\nabla\tau(x)\rvert^2=\frac{1}{c^2(x)}," title="\lvert\nabla\tau(x)\rvert^2=\frac{1}{c^2(x)}," /></p>  
and the factored eikonal equation  
<p align="center"><img src="https://latex.codecogs.com/svg.latex?\lvert(\tau_0\nabla\tau_1&plus;\tau_1\nabla\tau_0)(x)\rvert^2=\frac{1}{c^2(x)}," title="\lvert(\tau_0\nabla\tau_1+\tau_1\nabla\tau_0)(x)\rvert^2=\frac{1}{c^2(x)}," /></p>  
where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\tau_0(x)=\lvert&space;x-x_s\rvert." title="\tau_0(x)=\lvert x-x_s\rvert" />

## References
- J. Sethian. Fast marching methods. SIAM Review, 41(2):199-235, 1999. doi: 10.1137/S0036144598347059. URL https://doi.org/10.1137/S0036144598347059
- Eran Treister and Eldad Haber. A fast marching algorithm for the factored eikonal equation. Journal of Computational Physics, 324:210-225, 2016.


## Requirements

- Python 3
- numpy version 1.7 or higher
- C++11 compiler


## Installation

### Installation from PyPi:

```bash
pip install eikonalfm
```

### Manual install from the repository:

```bash
git clone https://github.com/kevinganster/eikonalfm.git
cd eikonalfm
pip install .
```

or

```bash
pip install git+https://github.com/kevinganster/eikonalfm.git
```


## Examples

```python
import numpy as np
import eikonalfm

c = np.ones((100, 100))
x_s = (0, 0)
dx = (1.0, 1.0)
order = 2

tau_fm = eikonalfm.fast_marching(c, x_s, dx, order)
tau1_ffm = eikonalfm.factored_fast_marching(c, x_s, dx, order)
```

Note that the source position `x_s` describes an index-vector.

To visualize the results, matplotlib (https://pypi.org/project/matplotlib/) can be used, for example:

```python
import matplotlib.pyplot as plt

# for the distance-function 'x_s' also describes an index-vector
tau0 = eikonalfm.distance(tau1_ffm.shape, dx, x_s, indexing="ij")
plt.contourf(tau0 * tau1_ffm)
plt.show()
```
