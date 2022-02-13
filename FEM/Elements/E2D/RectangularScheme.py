"""Define the rectangular scheme used by rectangular elements
"""


import numpy as np


class RectangularScheme():
    """Generates a rectangular integration scheme

    Args:
        n (int): Number of gauss points
    """

    def __init__(self, n: int, **kargs) -> None:
        """Generates a rectangular integration scheme

        Args:
            n (int): Number of gauss points
        """
        self.Tj = [
            lambda s: np.array([s, -1*(s-s+1)]),
            lambda s: np.array([1*(s-s+1), s]),
            lambda s: np.array([-1*s, 1*(s-s+1)]),
            lambda s: np.array([-1*(s-s+1), -1*s]),
        ]
        _Z, _W = np.polynomial.legendre.leggauss(n)
        self.Z = []
        self.W = []
        for i, z in enumerate(_Z):
            for j, n in enumerate(_Z):
                self.Z.append([z, n])
                self.W.append(_W[i]*_W[j])

        self.Z = np.array(self.Z)
        self.W = np.array(self.W)
        self.center = np.array([[0.0, 0.0]])
        self.domain = np.array(
            [[-1, -1], [1, -1], [1, 1], [-1, 1]] + self.Z.tolist())
