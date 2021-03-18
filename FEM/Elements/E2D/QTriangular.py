from .Element2D import *
from .TriangularScheme import *


class QTriangular(Element2D, TriangularScheme):
    def __init__(self, coords, gdl, n=3):
        coords = np.array(coords)
        _coords = np.array([coords[i] for i in range(3)])
        Element2D.__init__(self, coords, _coords, gdl)
        TriangularScheme.__init__(self, n)

    def psis(self, z):
        return np.array([
            2.0*(z[0]+z[1]-1.0)*(z[0]+z[1]-0.5),
            2.0*z[0]*(z[0]-0.5),
            2.0*z[1]*(z[1]-0.5),
            -4.0*(z[0]+z[1]-1.0)*(z[0]),
            4.0*z[0]*z[1],
            -4.0*z[1]*(z[0]+z[1]-1.0)]).T

    def dpsis(self, z):
        return np.array([
            [4.0*z[0]+4.0*z[1]-3.0, 4.0*z[1]+4.0*z[0]-3.0],
            [4.0*z[0]-1.0, 0*z[0]],
            [0*z[0], 4.0*z[1]-1.0],
            [-8.0*z[0]-4.0*(z[1]-1.0), -4.0*z[0]],
            [4.0*z[1], 4.0*z[0]],
            [-4.0*z[1], -8.0*z[1]-4.0*z[0]+4.0]
        ])
