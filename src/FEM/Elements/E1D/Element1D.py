"""1D Elements general class
"""


from ..Element import Element
from .LinearScheme import LinearScheme
import numpy as np
import matplotlib.pyplot as plt


class Element1D(Element, LinearScheme):
    """Create a 1D Element

    Args:
        coords (np.ndarray): Coordinates matrix
        gdl (np.ndarray): Degree of freedom matrix
        n (int): Number of Gauss Points used in integration
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int, **kargs) -> None:
        """Create a 1D Element

        Args:
            coords (np.ndarray): Coordinates matrix
            gdl (np.ndarray): Degree of freedom matrix
            n (int): Number of Gauss Points used in integration
        """

        coords = np.array(coords).reshape([len(coords), 1])
        _coords = np.array([coords[0][0], coords[-1][0]])
        LinearScheme.__init__(self, n)
        Element.__init__(self, coords, _coords, gdl, **kargs)

    def draw(self) -> None:
        """Create a element graph over domain

        """
        _z = self.domain
        _x, _p = self.T(_z)
        _y = 0
        for i in range(self.n):
            plt.plot(_x, _p[:, i], '--', label=r'$\psi_{'+format(i)+r'}$')
            _y += _p[:, i]
        plt.plot(_x, [0]*len(_x), '-', color='black', label='Element')
        plt.plot(self.coords.T[0], [0]*len(self.coords.T[0]),
                 'o', color='blue', label='Nodes')
        plt.plot(_x, _y, '-.', label=r'$\sum_{i=0}^{n}\psi$')
        plt.legend()
        plt.grid()

    def jacobianGraph(self) -> None:
        """Jacobian is constant in lineal elements
        """
        pass

    def isInside(self, x: np.ndarray) -> np.ndarray:
        """Test if a given points is inside element domain

        Args:
            x (np.ndarray): Point to be tested

        Returns:
            np.ndarray: Bolean array of test result
        """
        return (x >= self.coords.T[0][0])*(x <= self.coords.T[0][-1])
