# Eikonal Fast Marching

eikonalfm is a Python (C++-)extension which implements the fast marching method for the factored (and basic, i.e non-factored) eikonal equation  
```| grad tau(x) |^2 = 1 / c(x)^2.```


# Installation

Installation from PyPi:(TODO)  
```
pip install factoredeikonal
```

Manual install from the repository:  
```
git clone https://github.com/Daarknes/factoredeikonal.git
cd factoredeikonal
python setup.py
```


# Examples

```
import numpy as np
import eikonalfm

c = np.ones((100, 100))
x0 = np.array([0, 0])
dx = np.array([1.0, 1.0])
order = 2

tau_fm = eikonalfm.fast_marching(c, x0, dx, order)
tau_ffm = eikonalfm.factored_fast_marching(c, x0, dx, order)
```
