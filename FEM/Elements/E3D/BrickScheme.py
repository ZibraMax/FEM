"""Define the brick scheme used by brick elements
"""


import numpy as np


class BrickScheme():
    """Generate a brick integration scheme

    Args:
        n (int): Number of gauss points
    """

    def __init__(self, n: int, **kargs) -> None:
        """Generate a brick integration scheme

        Args:
            n (int): Number of gauss points
        """

        _Z, _W = np.polynomial.legendre.leggauss(n)
        self.Z = []
        self.W = []
        for i, z in enumerate(_Z):
            for j, n in enumerate(_Z):
                for k, g in enumerate(_Z):
                    self.Z.append([z, n, g])
                    self.W.append(_W[i]*_W[j]*_W[k])

        self.Z = np.array(self.Z)
        self.W = np.array(self.W)
        self.center = np.array([[0.0, 0.0, 0.0]])
        self.domain = np.array(
            [[-1, -1, -1],
             [1, -1, -1],
             [1, 1, -1],
             [-1, 1, -1],
             [-1, -1, 1],
             [1, -1, 1],
             [1, 1, 1],
             [-1, 1, 1]] + self.Z.tolist())
