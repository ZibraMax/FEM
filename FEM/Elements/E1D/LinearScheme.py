"""Defines a linear scheme for lineal elements
"""

import numpy as np


class LinearScheme():
    """Creates a Linear Scheme

    Args:
        n (int): Number of gauss points
    """

    def __init__(self, n: int) -> None:
        """Creates a Linear Scheme

        Args:
            n (int): Number of gauss points
        """

        self.Z, self.W = np.polynomial.legendre.leggauss(n)
        self.domain = np.array([[-1] + self.Z.tolist() + [1]])[0]
