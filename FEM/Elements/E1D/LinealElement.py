from .Element1D import *


class LinealElement(Element1D):
    def __init__(self, coords, gdl, n=3):
        Element1D.__init__(self, coords, gdl, n)

    def psis(self, z):
        return np.array([0.5*(1.0-z), 0.5*(1.0+z)]).T

    def dpsis(self, z):
        kernel = 1+(z-z)
        return np.array([[-0.5*kernel], [0.5*kernel]])
