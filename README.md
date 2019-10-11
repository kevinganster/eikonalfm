# Eikonal Fast Marching

eikonalfm is a Python (C++) extension which implements the Fast Marching method for the eikonal equation  
```| grad(tau)(x) |^2 = 1 / c^2(x)```,  
and the factored eikonal equation  
```| (tau0 grad(tau1) + tau1 grad(tau0))(x) |^2 = 1 / c^2(x)```,  
where `tau0(x) = | x - x_s |`.

## References
- J. Sethian. Fast marching methods. SIAM Review, 41(2):199-235, 1999. doi: 10.1137/S0036144598347059. URL https://doi.org/10.1137/S0036144598347059
- Eran Treister and Eldad Haber. A fast marching algorithm for the factored eikonal equation. Journal of Computational Physics, 324:210-225, 2016.


# Requirements

- Python 3
- numpy version 1.7 or above


# Installation

Installation from PyPi:  
```
pip install eikonalfm
```

Manual install from the repository:  
```
git clone https://github.com/Daarknes/eikonalfm.git
cd eikonalfm
python setup.py
```


# Examples

```
import numpy as np
import eikonalfm

c = np.ones((100, 100))
x_s = np.array([0, 0])
dx = np.array([1.0, 1.0])
order = 2

tau_fm = eikonalfm.fast_marching(c, x_s, dx, order)
tau_ffm = eikonalfm.factored_fast_marching(c, x_s, dx, order)
```

Note that the source position `x_s` describes an index-vector.
