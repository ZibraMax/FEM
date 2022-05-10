"""Defines a general 2D element
"""


from ..Element import Element
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpltPath


class Element2D(Element):

    """Create a 2D element

    Args:
        coords (np.ndarray): Element coordinate matrix
        _coords (np.ndarray): Element coordinate matrix for graphical interface purposes
        gdl (np.ndarray): Degree of freedom matrix
    """

    def __init__(self, coords: np.ndarray, _coords: np.ndarray, gdl: np.ndarray, **kargs) -> None:
        """Create a 2D element

        Args:
            coords (np.ndarray): Element coordinate matrix
            _coords (np.ndarray): Element coordinate matrix for graphical interface purposes
            gdl (np.ndarray): Degree of freedom matrix
        """

        Element.__init__(self, coords, _coords, gdl, **kargs)
        self._coordsg = np.array(
            self._coords.tolist()+[self._coords[0].tolist()])
        for i, e in enumerate(self.borders):
            delta = self._coordsg[i+1]-self._coordsg[i]
            delta[0] *= -1
            delta = delta[::-1]
            delta = delta/np.linalg.norm(delta)
            e.nx = delta[0]
            e.ny = delta[1]

    def draw(self) -> None:
        """Create a graph of element"""

        _z = self.domain
        _x, _p = self.T(_z.T)
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        l = []
        l.append('Element')
        l.append('Nodes')
        for i in range(self.n):
            surf = ax.plot_trisurf(*_x.T, _p[:, i], alpha=0.3)
            surf._facecolors2d = surf._facecolor3d
            surf._edgecolors2d = surf._edgecolor3d
            l.append(r'$\psi_{'+format(i)+r'}$')
        __coords = np.array(self._coords.tolist()+[self._coords[0].tolist()]).T
        ax.plot(*__coords, [0]*len(__coords.T), '-', color='black')
        ax.plot(*self.coords.T, [0]*len(self.coords), 'o', color='blue')
        ax.legend(l)

    def jacobianGraph(self) -> None:
        """Create the determinant jacobian graph
        """

        _z = self.domain
        _x, _p = self.T(_z.T)
        _j = self.J(_z.T)[0]
        __j = np.linalg.det(_j)
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        l = []
        surf = ax.plot_trisurf(*_x.T, __j, cmap='magma')
        surf._facecolors2d = surf._facecolor3d
        surf._edgecolors2d = surf._edgecolor3d
        l.append('Element')
        l.append('Nodes')
        l.append(r'$|J|$')
        fig.colorbar(surf)
        __coords = np.array(self._coords.tolist()+[self._coords[0].tolist()]).T
        ax.plot(*__coords, [0]*len(__coords.T), '-', color='black')
        ax.plot(*self.coords.T, [0]*len(self.coords), 'o', color='blue')
        ax.legend(l)

    def isInside(self, x: np.ndarray) -> np.ndarray:
        """Test if a given points is inside element domain

        Args:
            x (np.ndarray): Point to be tested

        Returns:
            np.ndarray: Bolean array of test result
        """
        path = mpltPath.Path(self._coords[:, :2])
        inside2 = path.contains_points([x])
        return inside2[0]
