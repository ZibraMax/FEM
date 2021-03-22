"""Define the rectangular scheme used by rectangular elements
"""


import numpy as np


class RectangularScheme():
    """Generates a rectangular integration scheme

    Args:
        n (int): Number of gauss points
    """

    def __init__(self, n: int) -> None:
        """Generates a rectangular integration scheme

        Args:
            n (int): Number of gauss points
        """

        DENSIDAD = 10
        XX, YY = np.meshgrid(np.linspace(-1, 1, DENSIDAD),
                             np.linspace(-1, 1, DENSIDAD))
        _z = XX.reshape([DENSIDAD**2, 1])[:, 0]
        _n = YY.reshape([DENSIDAD**2, 1])[:, 0]
        _Z, _W = np.polynomial.legendre.leggauss(n)
        self.Z = []
        self.W = []
        for i, z in enumerate(_Z):
            for j, n in enumerate(_Z):
                self.Z.append([z, n])
                self.W.append(_W[i]*_W[j])

        self.Z = np.array(self.Z)
        self.W = np.array(self.W)
        self.domain = np.array(
            [[-1, -1], [1, -1], [1, 1], [-1, 1]] + self.Z.tolist())
