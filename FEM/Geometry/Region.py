import numpy as np

from FEM import Utils
from ..Elements import Quadrilateral
from ..Utils import isBetween


class Region():
    """docstring for Region
    """
    # TODO hacer los docstrigns

    def __init__(self, coords):
        if not isinstance(coords, np.ndarray):
            coords = np.array(coords)
        self.coords = coords
        self.nodes = []

    def setNodesOfRegion(self, geometry, tol=10**(-5)):
        self.nodes = []
        for i, p in enumerate(geometry.gdls):
            if self.isBetween(p, tol):
                self.nodes.append(i)
        self.nodes = np.array(self.nodes)


class Region1D(Region):
    """docstring for Region1D
    """

    def __init__(self, coords):
        Region.__init__(self, coords)

    def isBetween(self, p, tol: float = 1*10**(-5)):
        return isBetween(self.coords[0], self.coords[1], p, tol)


class Region2D(Region):
    """docstring for Region2D
    """

    def __init__(self, coords):
        self.e = Quadrilateral(coords, np.array([[-1, -1, -1, -1]]), n=1)
        self.center, _ = self.e.T(self.e.center.T)
        _j, _ = self.e.J(self.e.center.T)
        self.n = np.cross(_j[:, 0, :].T, _j[:, 1, :].T, axis=0)  # <A,B,C>
        self.nnorm = np.linalg.norm(self.n)
        # Ax + By + Cz + D = 0
        self.D = -np.dot(self.n.flatten(), self.center.flatten())
        Region.__init__(self, coords)

    def pointToPlaneDistance(self, p):
        num = abs(np.dot(self.n.T.flatten(), p.flatten())+self.D)
        return num/self.nnorm

    def isBetween(self, p, tol: float = 1*10**(-5)):
        d = self.pointToPlaneDistance(p)
        return d <= tol
