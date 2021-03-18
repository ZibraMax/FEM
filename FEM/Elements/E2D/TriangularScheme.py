import numpy as np


class TriangularScheme():
    def __init__(self, n):
        # self.Z,self.W = np.polynomial.legendre.leggauss(n)
        # w = np.outer((1 + points) * weights, weights) / 4
        # x = np.outer(1 - points, np.ones(n)) / 2
        # y = np.outer(1 + points, 1 - points) / 4
        # points = np.array([x.flatten(), y.flatten()])
        # weights = w.flatten()
        # points = np.array([points[0], points[1], 1 - points[0] - points[1]])
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
        self.domain = np.array([[0, 0], [1, 0], [0, 1]] + self.Z.tolist())
