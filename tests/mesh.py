import numpy as np
import sys


class Mesh2D:
    def __init__(self, x1_range, x2_range, nx1=None, nx2=None):
        if nx1 is None and nx2 is None:
            raise ValueError

        if nx2 is None:
            x1, dx1 = np.linspace(x1_range[0], x1_range[1], num=nx1, retstep=True)
            nx2 = int((x2_range[1] - x2_range[0]) / dx1) + 1
            x2, dx2 = np.linspace(x2_range[0], x2_range[1], num=nx2, retstep=True)
        else:
            x2, dx2 = np.linspace(x2_range[0], x2_range[1], num=nx2, retstep=True)
            nx1 = int((x1_range[1] - x1_range[0]) / dx2) + 1
            x1, dx1 = np.linspace(x1_range[0], x1_range[1], num=nx1, retstep=True)

        self.n = (nx1, nx2)
        self.dx = np.array([dx1, dx2])
        self.X = np.meshgrid(x1, x2, indexing="ij")
        self.origin = np.array([self.X[0][0, 0], self.X[1][0, 0]])

    def index_to_pos(self, x):
        return self.origin + x * self.dx

    def pos_to_index(self, x):
        return ((x - self.origin + sys.float_info.epsilon) / self.dx).astype(np.int)

    def __str__(self):
        endpoint = (self.X[0][-1, -1], self.X[1][-1, -1])
        return """Mesh2D:
    area: [{o[0]}, {e[0]}] x [{o[1]}, {e[1]}]
    n = {n}
    dx = ({dx[0]}, {dx[1]})""".format(o=self.origin, e=endpoint, n=self.n, dx=self.dx)