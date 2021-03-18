import numpy as np


class LinearScheme():
    def __init__(self, n):
        self.Z, self.W = np.polynomial.legendre.leggauss(n)
        self.domain = np.array([[-1] + self.Z.tolist() + [1]])[0]

    def integrate(self, f):
        integral = 0
        for w, z in zip(self.W, self.Z):
            integral += f(z)*w
        return integral
