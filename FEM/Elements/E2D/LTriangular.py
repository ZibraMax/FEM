from .Element2D import *
from .TriangularScheme import *


class LTriangular(Element2D, TriangularScheme):
    def __init__(self, coords, gdl, n=2):
        coords = np.array(coords)
        Element2D.__init__(self, coords, coords, gdl)
        TriangularScheme.__init__(self, n)

    def psis(self, z):
        return np.array([
            1.0-z[0]-z[1],
            z[0],
            z[1]]).T

    def dpsis(self, z):
        kernell = (z[0]-z[0])
        return np.array([
            [-1.0*(1+kernell), -1.0*(1+kernell)],
            [1.0*(1+kernell), 0.0*(1+kernell)],
            [0.0*(1+kernell), 1.0*(1+kernell)]
        ])
