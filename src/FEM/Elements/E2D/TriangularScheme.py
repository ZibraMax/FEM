"""Define the triangular scheme used by triangular elements
"""


import numpy as np


class TriangularScheme():
    """Generate a triangular integration scheme

    Args:
        n (int): Number of gauss points
    """

    def __init__(self, n: int, **kargs) -> None:
        """Generate a triangular integration scheme

        Args:
            n (int): Number of gauss points
        """

        # self.Z,self.W = np.polynomial.legendre.leggauss(n)
        # w = np.outer((1 + points) * weights, weights) / 4
        # x = np.outer(1 - points, np.ones(n)) / 2
        # y = np.outer(1 + points, 1 - points) / 4
        # points = np.array([x.flatten(), y.flatten()])
        # weights = w.flatten()
        # points = np.array([points[0], points[1], 1 - points[0] - points[1]])

        self.Tj = [
            lambda s: np.array([(0.5*s+0.5), (s-s)]),
            lambda s: np.array([1-(0.5*s+0.5), (0.5*s+0.5)]),
            lambda s: np.array([(s-s), 1-(0.5*s+0.5)]),
        ]

        A0 = 1/3
        A1 = 0.059715871789770
        A2 = 0.797426985353087
        B1 = 0.470142064105115
        B2 = 0.101286507323456
        W0 = 0.1125
        W1 = 0.066197076394253
        W2 = 0.062969590272413
        X = [A0, A1, B1, B1, B2, B2, A2]
        Y = [A0, B1, A1, B1, A2, B2, B2]
        W = [W0, W1, W1, W1, W2, W2, W2]
        self.Z = np.array([X, Y]).T
        self.W = np.array(W)
        self.center = np.array([[1.0, 1.0]])/3.0
        self.domain = np.array([[0, 0], [1, 0], [0, 1]] + self.Z.tolist())
