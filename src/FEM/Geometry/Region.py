"""Defines geometry region. Regions are used to set border conditions and loads to 2D and 3D problems"""

import numpy as np

from ..Elements import Quadrilateral, Serendipity, QTriangular, LTriangular
from ..Utils import isBetween


class Region():
    """Creates a general region

    Args:
        coords (np.ndarray): Coordinates matrix of therefion. Each row is a coordinate, each column is a dimension.
    """

    def __init__(self, coords: np.ndarray, desc: str = '') -> None:
        """Creates a general region

        Args:
            coords (np.ndarray): Coordinates matrix of therefion. Each row is a coordinate, each column is a dimension.
            desc (str, optional): Region description.
        """
        if not isinstance(coords, np.ndarray):
            coords = np.array(coords)
        self.coords = coords
        self.nodes = []
        self.description = desc

    def setNodesOfRegion(self, geometry: 'Geometry', tol: float = 10**(-5)) -> None:
        """Calculates the nodes of the geometry which are inside the region.

        Args:
            geometry (Geometry): Input geometry
            tol (float, optional): Near tolerance. Defaults to 10**(-5).
        """
        self.nodes = []
        for i, p in enumerate(geometry.gdls):
            if self.isBetween(p, tol):
                self.nodes.append(i)
        self.nodes = np.array(self.nodes)


class Region1D(Region):
    """Creates a line region (1D element)

    Args:
        coords (np.ndarray): Coordinate matrix. Must be of two rows and 2 or 3 columns. 2 columns for 2D region, 3 columns for 3D region.
    """

    def __init__(self, coords: np.ndarray, **kargs) -> None:
        """Creates a line region (1D element)

        Args:
            coords (np.ndarray): Coordinate matrix. Must be of two rows and 2 or 3 columns. 2 columns for 2D region, 3 columns for 3D region.
        """
        Region.__init__(self, coords, **kargs)

    def isBetween(self, p: np.ndarray, tol: float = 1*10**(-5)) -> bool:
        """Check if a given point is inside the region.

        Args:
            p (np.ndarray): Point to be tested
            tol (float, optional): Tolerance for check. Defaults to 1*10**(-5).

        Returns:
            bool: True if the point is inside the region
        """
        return isBetween(self.coords[0], self.coords[1], p, tol)


class Region2D(Region):
    """Creates a square region (2D element)

    Args:
        coords (np.ndarray): Coordinate matrix. Must be of four rows and 3 columns.
    """

    def __init__(self, coords: np.ndarray, **kargs) -> None:
        """Creates a 2D region using a 2D Element. The number of coordinates definesthe type of element.

        Args:
            coords (np.ndarray): Coordinate matrix. Must be of four rows and 3 columns.
        """
        ndim = len(coords[0])
        lc = len(coords)
        if ndim == 2:
            coordsn = np.zeros([lc, 3])
            coordsn[:, :-1] = coords
            coords = coordsn
        if lc == 3:
            ELE = LTriangular
        elif lc == 4:
            ELE = Quadrilateral
        elif lc == 6:
            ELE = QTriangular
        elif lc == 8:
            ELE = Serendipity
        self.e = ELE(coords, np.array([[-1]*lc]), n=1, fast=True, border=True)
        self.center, _ = self.e.T(self.e.center.T)
        _j, _ = self.e.J(self.e.center.T)
        self.n = np.cross(_j[:, 0].T, _j[:, 1].T, axis=0)
        self.nnorm = np.linalg.norm(self.n)
        # Ax + By + Cz + D = 0
        self.D = -np.dot(self.n.flatten(), self.center.flatten())
        Region.__init__(self, coords, **kargs)

    def pointToPlaneDistance(self, p: np.ndarray) -> float:
        """Calculates the distance from a given point to the region.

        Args:
            p (np.ndarray): Point to be tested

        Returns:
            float: Distance between the plane and the point
        """
        num = abs(np.dot(self.n.T.flatten(), p.flatten())+self.D)
        return num/self.nnorm

    def isBetween(self, p: np.ndarray, tol: float = 1*10**(-5)) -> bool:
        """Check if a given point is inside the region.

        Args:
            p (np.ndarray): Point to be tested
            tol (float, optional): Tolerance for check. Defaults to 1*10**(-5).

        Returns:
            bool: True if the point is inside the region
        """
        ndim = len(p)
        if ndim == 2:
            return self.e.isInside(p)
        d = self.pointToPlaneDistance(p)
        return d <= tol
