import numpy as np
import time
# import peikonal
import eikonalfm
np.set_printoptions(precision=2, linewidth=1e6)


#np.random.seed = 42
c = 1 * np.ones((1000, 1000))
# c = np.array([
#     [6, 1, 1],
#     [6, 1, 6],
#     [6, 6, 6]
# ])
x0 = np.array([500, 500])
dx = np.array([0.1, 0.1])
order = 2


start = time.perf_counter()
tau = eikonalfm.fast_marching(c, x0, dx, order)
print("own version took {:.8f} s".format(time.perf_counter() - start))
print(np.array2string(tau, separator="    "))
 
 
print("-"*50)
start = time.perf_counter()
tau = eikonalfm.factored_fast_marching(c, x0, dx, order)
print("factored version took {:.8f} s".format(time.perf_counter() - start))
print(np.array2string(tau, separator="    "))


print("-"*50)
# comparison to skfmm
import skfmm
phi = np.ones_like(c)
phi[tuple(x0)] = 0
    
start = time.perf_counter()
tau2 = skfmm.travel_time(phi, c, dx=dx, order=order)
print("skfmm version took {:.8f} s".format(time.perf_counter() - start))
print(np.array2string(tau2, separator="    "))

# print(tau1 - tau2)