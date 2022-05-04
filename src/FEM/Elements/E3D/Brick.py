"""
*********************************************
BRICK ELEMENTS
*********************************************

Defines the lagrangian first order and second order brick elements

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/docs/source/Bricks.png
   :align: center

   Brick (8 nodes) and BrickO2 (20 nodes). (Reddy, 2005)

Brick
#####

First order 8 node brick element

Shape Functions
---------------

.. math::
    \\Psi_i=\\frac{1}{8}\\left[\\begin{matrix}\\left(1 - x\\right) \\left(1 - y\\right) \\left(1 - z\\right)\\\\\\left(1 - y\\right) \\left(1 - z\\right) \\left(x + 1\\right)\\\\\\left(1 - z\\right) \\left(x + 1\\right) \\left(y + 1\\right)\\\\\\left(1 - x\\right) \\left(1 - z\\right) \\left(y + 1\\right)\\\\\\left(1 - x\\right) \\left(1 - y\\right) \\left(z + 1\\right)\\\\\\left(1 - y\\right) \\left(x + 1\\right) \\left(z + 1\\right)\\\\\\left(x + 1\\right) \\left(y + 1\\right) \\left(z + 1\\right)\\\\\\left(1 - x\\right) \\left(y + 1\\right) \\left(z + 1\\right)\\end{matrix}\\right]

Shape Functions Derivatives
---------------------------

.. math::
    \\frac{\\partial \\Psi_i}{\\partial x_j}=\\frac{1}{8}\\left[\\begin{matrix}\\left(1 - z\\right) \\left(y - 1\\right) & \\left(1 - z\\right) \\left(x - 1\\right) & - \\left(1 - x\\right) \\left(1 - y\\right)\\\\\\left(1 - y\\right) \\left(1 - z\\right) & \\left(1 - z\\right) \\left(- x - 1\\right) & - \\left(1 - y\\right) \\left(x + 1\\right)\\\\\\left(1 - z\\right) \\left(y + 1\\right) & \\left(1 - z\\right) \\left(x + 1\\right) & - \\left(x + 1\\right) \\left(y + 1\\right)\\\\\\left(1 - z\\right) \\left(- y - 1\\right) & \\left(1 - x\\right) \\left(1 - z\\right) & - \\left(1 - x\\right) \\left(y + 1\\right)\\\\\\left(1 - y\\right) \\left(- z - 1\\right) & - \\left(1 - x\\right) \\left(z + 1\\right) & \\left(1 - x\\right) \\left(1 - y\\right)\\\\\\left(1 - y\\right) \\left(z + 1\\right) & - \\left(x + 1\\right) \\left(z + 1\\right) & \\left(1 - y\\right) \\left(x + 1\\right)\\\\\\left(y + 1\\right) \\left(z + 1\\right) & \\left(x + 1\\right) \\left(z + 1\\right) & \\left(x + 1\\right) \\left(y + 1\\right)\\\\- \\left(y + 1\\right) \\left(z + 1\\right) & \\left(1 - x\\right) \\left(z + 1\\right) & \\left(1 - x\\right) \\left(y + 1\\right)\\end{matrix}\\right]

BrickO2
#######

Second order 20 node brick element

Shape Functions
---------------

.. math::
    \\Psi_i = \\frac{1}{8}\\left[\\begin{matrix}\\left(1 - x\\right) \\left(1 - y\\right) \\left(1 - z\\right) \\left(- x - y - z - 2\\right)\\\\\\left(1 - y\\right) \\left(1 - z\\right) \\left(x + 1\\right) \\left(x - y - z - 2\\right)\\\\\\left(1 - z\\right) \\left(x + 1\\right) \\left(y + 1\\right) \\left(x + y - z - 2\\right)\\\\\\left(1 - x\\right) \\left(1 - z\\right) \\left(y + 1\\right) \\left(- x + y - z - 2\\right)\\\\\\left(1 - x\\right) \\left(1 - y\\right) \\left(z + 1\\right) \\left(- x - y + z - 2\\right)\\\\\\left(1 - y\\right) \\left(x + 1\\right) \\left(z + 1\\right) \\left(x - y + z - 2\\right)\\\\\\left(x + 1\\right) \\left(y + 1\\right) \\left(z + 1\\right) \\left(x + y + z - 2\\right)\\\\\\left(1 - x\\right) \\left(y + 1\\right) \\left(z + 1\\right) \\left(- x + y + z - 2\\right)\\\\\\left(1 - y\\right) \\left(1 - z\\right) \\left(2 - 2 x^{2}\\right)\\\\\\left(1 - y^{2}\\right) \\left(1 - z\\right) \\left(2 x + 2\\right)\\\\\\left(1 - z\\right) \\left(2 - 2 x^{2}\\right) \\left(y + 1\\right)\\\\\\left(1 - y^{2}\\right) \\left(1 - z\\right) \\left(2 - 2 x\\right)\\\\\\left(1 - y\\right) \\left(1 - z^{2}\\right) \\left(2 - 2 x\\right)\\\\\\left(1 - y\\right) \\left(1 - z^{2}\\right) \\left(2 x + 2\\right)\\\\\\left(1 - z^{2}\\right) \\left(2 x + 2\\right) \\left(y + 1\\right)\\\\\\left(1 - z^{2}\\right) \\left(2 - 2 x\\right) \\left(y + 1\\right)\\\\\\left(1 - y\\right) \\left(2 - 2 x^{2}\\right) \\left(z + 1\\right)\\\\\\left(1 - y^{2}\\right) \\left(2 x + 2\\right) \\left(z + 1\\right)\\\\\\left(2 - 2 x^{2}\\right) \\left(y + 1\\right) \\left(z + 1\\right)\\\\\\left(1 - y^{2}\\right) \\left(2 - 2 x\\right) \\left(z + 1\\right)\\end{matrix}\\right]

Shape Functions Derivatives
---------------------------

.. math::
    \\frac{\\partial \\Psi_i}{\\partial x_j} = \\frac{1}{8}\\left[\\begin{matrix}- \\left(1 - x\\right) \\left(1 - y\\right) \\left(1 - z\\right) + \\left(1 - z\\right) \\left(y - 1\\right) \\left(- x - y - z - 2\\right) & - \\left(1 - x\\right) \\left(1 - y\\right) \\left(1 - z\\right) + \\left(1 - z\\right) \\left(x - 1\\right) \\left(- x - y - z - 2\\right) & - \\left(1 - x\\right) \\left(1 - y\\right) \\left(1 - z\\right) - \\left(1 - x\\right) \\left(1 - y\\right) \\left(- x - y - z - 2\\right)\\\\\\left(1 - y\\right) \\left(1 - z\\right) \\left(x + 1\\right) + \\left(1 - y\\right) \\left(1 - z\\right) \\left(x - y - z - 2\\right) & - \\left(1 - y\\right) \\left(1 - z\\right) \\left(x + 1\\right) + \\left(1 - z\\right) \\left(- x - 1\\right) \\left(x - y - z - 2\\right) & - \\left(1 - y\\right) \\left(1 - z\\right) \\left(x + 1\\right) - \\left(1 - y\\right) \\left(x + 1\\right) \\left(x - y - z - 2\\right)\\\\\\left(1 - z\\right) \\left(x + 1\\right) \\left(y + 1\\right) + \\left(1 - z\\right) \\left(y + 1\\right) \\left(x + y - z - 2\\right) & \\left(1 - z\\right) \\left(x + 1\\right) \\left(y + 1\\right) + \\left(1 - z\\right) \\left(x + 1\\right) \\left(x + y - z - 2\\right) & - \\left(1 - z\\right) \\left(x + 1\\right) \\left(y + 1\\right) - \\left(x + 1\\right) \\left(y + 1\\right) \\left(x + y - z - 2\\right)\\\\- \\left(1 - x\\right) \\left(1 - z\\right) \\left(y + 1\\right) + \\left(1 - z\\right) \\left(- y - 1\\right) \\left(- x + y - z - 2\\right) & \\left(1 - x\\right) \\left(1 - z\\right) \\left(y + 1\\right) + \\left(1 - x\\right) \\left(1 - z\\right) \\left(- x + y - z - 2\\right) & - \\left(1 - x\\right) \\left(1 - z\\right) \\left(y + 1\\right) - \\left(1 - x\\right) \\left(y + 1\\right) \\left(- x + y - z - 2\\right)\\\\- \\left(1 - x\\right) \\left(1 - y\\right) \\left(z + 1\\right) + \\left(1 - y\\right) \\left(- z - 1\\right) \\left(- x - y + z - 2\\right) & - \\left(1 - x\\right) \\left(1 - y\\right) \\left(z + 1\\right) - \\left(1 - x\\right) \\left(z + 1\\right) \\left(- x - y + z - 2\\right) & \\left(1 - x\\right) \\left(1 - y\\right) \\left(z + 1\\right) + \\left(1 - x\\right) \\left(1 - y\\right) \\left(- x - y + z - 2\\right)\\\\\\left(1 - y\\right) \\left(x + 1\\right) \\left(z + 1\\right) + \\left(1 - y\\right) \\left(z + 1\\right) \\left(x - y + z - 2\\right) & - \\left(1 - y\\right) \\left(x + 1\\right) \\left(z + 1\\right) - \\left(x + 1\\right) \\left(z + 1\\right) \\left(x - y + z - 2\\right) & \\left(1 - y\\right) \\left(x + 1\\right) \\left(z + 1\\right) + \\left(1 - y\\right) \\left(x + 1\\right) \\left(x - y + z - 2\\right)\\\\\\left(x + 1\\right) \\left(y + 1\\right) \\left(z + 1\\right) + \\left(y + 1\\right) \\left(z + 1\\right) \\left(x + y + z - 2\\right) & \\left(x + 1\\right) \\left(y + 1\\right) \\left(z + 1\\right) + \\left(x + 1\\right) \\left(z + 1\\right) \\left(x + y + z - 2\\right) & \\left(x + 1\\right) \\left(y + 1\\right) \\left(z + 1\\right) + \\left(x + 1\\right) \\left(y + 1\\right) \\left(x + y + z - 2\\right)\\\\- \\left(1 - x\\right) \\left(y + 1\\right) \\left(z + 1\\right) - \\left(y + 1\\right) \\left(z + 1\\right) \\left(- x + y + z - 2\\right) & \\left(1 - x\\right) \\left(y + 1\\right) \\left(z + 1\\right) + \\left(1 - x\\right) \\left(z + 1\\right) \\left(- x + y + z - 2\\right) & \\left(1 - x\\right) \\left(y + 1\\right) \\left(z + 1\\right) + \\left(1 - x\\right) \\left(y + 1\\right) \\left(- x + y + z - 2\\right)\\\\- 4 x \\left(1 - y\\right) \\left(1 - z\\right) & \\left(2 - 2 x^{2}\\right) \\left(z - 1\\right) & \\left(2 - 2 x^{2}\\right) \\left(y - 1\\right)\\\\2 \\left(1 - y^{2}\\right) \\left(1 - z\\right) & - 2 y \\left(1 - z\\right) \\left(2 x + 2\\right) & \\left(2 x + 2\\right) \\left(y^{2} - 1\\right)\\\\- 4 x \\left(1 - z\\right) \\left(y + 1\\right) & \\left(1 - z\\right) \\left(2 - 2 x^{2}\\right) & \\left(2 - 2 x^{2}\\right) \\left(- y - 1\\right)\\\\- 2 \\left(1 - y^{2}\\right) \\left(1 - z\\right) & - 2 y \\left(1 - z\\right) \\left(2 - 2 x\\right) & \\left(2 - 2 x\\right) \\left(y^{2} - 1\\right)\\\\- 2 \\left(1 - y\\right) \\left(1 - z^{2}\\right) & \\left(2 - 2 x\\right) \\left(z^{2} - 1\\right) & - 2 z \\left(1 - y\\right) \\left(2 - 2 x\\right)\\\\2 \\left(1 - y\\right) \\left(1 - z^{2}\\right) & \\left(2 x + 2\\right) \\left(z^{2} - 1\\right) & - 2 z \\left(1 - y\\right) \\left(2 x + 2\\right)\\\\2 \\left(1 - z^{2}\\right) \\left(y + 1\\right) & \\left(1 - z^{2}\\right) \\left(2 x + 2\\right) & - 2 z \\left(2 x + 2\\right) \\left(y + 1\\right)\\\\- 2 \\left(1 - z^{2}\\right) \\left(y + 1\\right) & \\left(1 - z^{2}\\right) \\left(2 - 2 x\\right) & - 2 z \\left(2 - 2 x\\right) \\left(y + 1\\right)\\\\- 4 x \\left(1 - y\\right) \\left(z + 1\\right) & \\left(2 - 2 x^{2}\\right) \\left(- z - 1\\right) & \\left(1 - y\\right) \\left(2 - 2 x^{2}\\right)\\\\2 \\left(1 - y^{2}\\right) \\left(z + 1\\right) & - 2 y \\left(2 x + 2\\right) \\left(z + 1\\right) & \\left(1 - y^{2}\\right) \\left(2 x + 2\\right)\\\\- 4 x \\left(y + 1\\right) \\left(z + 1\\right) & \\left(2 - 2 x^{2}\\right) \\left(z + 1\\right) & \\left(2 - 2 x^{2}\\right) \\left(y + 1\\right)\\\\- 2 \\left(1 - y^{2}\\right) \\left(z + 1\\right) & - 2 y \\left(2 - 2 x\\right) \\left(z + 1\\right) & \\left(1 - y^{2}\\right) \\left(2 - 2 x\\right)\\end{matrix}\\right]

"""

from ..E2D.Quadrilateral import Quadrilateral
from ..E2D.Serendipity import Serendipity
from .Element3D import Element3D, np
from .BrickScheme import BrickScheme


class Brick(Element3D, BrickScheme):
    """Creates a 3D brick element

    Args:
        coords (np.ndarray): Node coordinates matrix
        gdl (np.ndarray): Degrees of freedom matrix
        n (int, optional): Number of gauss points used for integration. Defaults to 3.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 3, **kargs) -> None:
        """Creates a 3D brick element

        Args:
            coords (np.ndarray): Node coordinates matrix
            gdl (np.ndarray): Degrees of freedom matrix
            n (int, optional): Number of gauss points used for integration. Defaults to 3.
        """

        coords = np.array(coords)
        self.faces = [
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [4, 5, 6, 7],
            [3, 2, 1, 0],
            [2, 3, 7, 6],
            [4, 7, 3, 0]]
        self.face_element = Quadrilateral

        BrickScheme.__init__(self, n, **kargs)
        Element3D.__init__(self, coords, coords, gdl, **kargs)

    def psis(self, _z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """

        z = _z[0]
        n = _z[1]
        g = _z[2]
        return 1/8*np.array(
            [(1-z)*(1-n)*(1-g),
             (1+z)*(1-n)*(1-g),
             (1+z)*(1+n)*(1-g),
             (1-z)*(1+n)*(1-g),
             (1-z)*(1-n)*(1+g),
             (1+z)*(1-n)*(1+g),
             (1+z)*(1+n)*(1+g),
             (1-z)*(1+n)*(1+g)]).T

    def dpsis(self, _z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """
        x = _z[0]
        y = _z[1]
        z = _z[2]
        return 1/8*np.array(
            [[(y-1.0)*(1.0-z), (x-1)*(1-z), -(1-x)*(1-y)],
             [(1-y)*(1-z), (-1-x)*(1-z), -(1+x)*(1-y)],
             [(1+y)*(1-z), (1+x)*(1-z), -(1+x)*(1+y)],
             [(-1.0-y)*(1-z), (1-x)*(1-z), -(1-x)*(1+y)],
             [(1-y)*(-1-z), -(1-x)*(1+z), (1-x)*(1-y)],
             [(1-y)*(1+z), -(1+x)*(1+z), (1+x)*(1-y)],
             [(1+y)*(1+z), (1+x)*(1+z), (1+x)*(1+y)],
             [-(1+y)*(1+z), (1-x)*(1+z), (1-x)*(1+y)]])


class BrickO2(Element3D, BrickScheme):
    """Creates a 3D second order brick element.

    Args:
        coords (np.ndarray): Node coordinates matrix
        gdl (np.ndarray): Degrees of freedom matrix
        n (int, optional): Number of gauss points used for integration. Defaults to 3.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 3, **kargs) -> None:
        """Creates a 3D second order brick element.

        Args:
            coords (np.ndarray): Node coordinates matrix
            gdl (np.ndarray): Degrees of freedom matrix
            n (int, optional): Number of gauss points used for integration. Defaults to 3.
        """

        coords = np.array(coords)
        self.faces = [
            [0, 1, 5, 4, 8, 13, 16, 12],
            [1, 2, 6, 5, 9, 14, 17, 13],
            [4, 5, 6, 7, 16, 17, 18, 19],
            [3, 2, 1, 0, 10, 9, 8, 11],
            [2, 3, 7, 6, 10, 15, 18, 14],
            [4, 7, 3, 0, 19, 15, 11, 12]]
        self.face_element = Serendipity

        BrickScheme.__init__(self, n, **kargs)
        Element3D.__init__(self, coords, coords, gdl, **kargs)

    def psis(self, _z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """

        x = _z[0]
        y = _z[1]
        z = _z[2]
        return 1/8*np.array([
            (1-x)*(1-y)*(1-z)*(-x-y-z-2),
            (1+x)*(1-y)*(1-z)*(x-y-z-2),
            (1+x)*(1+y)*(1-z)*(x+y-z-2),
            (1-x)*(1+y)*(1-z)*(-x+y-z-2),
            (1-x)*(1-y)*(1+z)*(-x-y+z-2),
            (1+x)*(1-y)*(1+z)*(x-y+z-2),
            (1+x)*(1+y)*(1+z)*(x+y+z-2),
            (1-x)*(1+y)*(1+z)*(-x+y+z-2),
            2*(1-x**2)*(1-y)*(1-z),
            2*(1+x)*(1-y**2)*(1-z),
            2*(1-x**2)*(1+y)*(1-z),
            2*(1-x)*(1-y**2)*(1-z),
            2*(1-x)*(1-y)*(1-z**2),
            2*(1+x)*(1-y)*(1-z**2),
            2*(1+x)*(1+y)*(1-z**2),
            2*(1-x)*(1+y)*(1-z**2),
            2*(1-x**2)*(1-y)*(1+z),
            2*(1+x)*(1-y**2)*(1+z),
            2*(1-x**2)*(1+y)*(1+z),
            2*(1-x)*(1-y**2)*(1+z)
        ]).T

    def dpsis(self, _z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """
        x = _z[0]
        y = _z[1]
        z = _z[2]  # I JUST WANT THIS TO WORK PROPERLY
        return 1/8*np.array([[-(1 - x)*(1 - y)*(1 - z) + (1 - z)*(y - 1)*(-x - y - z - 2), -(1 - x)*(1 - y)*(1 - z) + (1 - z)*(x - 1)*(-x - y - z - 2), -(1 - x)*(1 - y)*(1 - z) - (1 - x)*(1 - y)*(-x - y - z - 2)],
                             [(1 - y)*(1 - z)*(x + 1) + (1 - y)*(1 - z)*(x - y - z - 2), -(1 - y)*(1 - z)*(x + 1) + (
                                 1 - z)*(-x - 1)*(x - y - z - 2), -(1 - y)*(1 - z)*(x + 1) - (1 - y)*(x + 1)*(x - y - z - 2)],
                             [(1 - z)*(x + 1)*(y + 1) + (1 - z)*(y + 1)*(x + y - z - 2), (1 - z)*(x + 1)*(y + 1) + (
                                 1 - z)*(x + 1)*(x + y - z - 2), -(1 - z)*(x + 1)*(y + 1) - (x + 1)*(y + 1)*(x + y - z - 2)],
                             [-(1 - x)*(1 - z)*(y + 1) + (1 - z)*(-y - 1)*(-x + y - z - 2), (1 - x)*(1 - z)*(y + 1) + (1 - x)
                              * (1 - z)*(-x + y - z - 2), -(1 - x)*(1 - z)*(y + 1) - (1 - x)*(y + 1)*(-x + y - z - 2)],
                             [-(1 - x)*(1 - y)*(z + 1) + (1 - y)*(-z - 1)*(-x - y + z - 2), -(1 - x)*(1 - y)*(z + 1) -
                              (1 - x)*(z + 1)*(-x - y + z - 2), (1 - x)*(1 - y)*(z + 1) + (1 - x)*(1 - y)*(-x - y + z - 2)],
                             [(1 - y)*(x + 1)*(z + 1) + (1 - y)*(z + 1)*(x - y + z - 2), -(1 - y)*(x + 1)*(z + 1) -
                              (x + 1)*(z + 1)*(x - y + z - 2), (1 - y)*(x + 1)*(z + 1) + (1 - y)*(x + 1)*(x - y + z - 2)],
                             [(x + 1)*(y + 1)*(z + 1) + (y + 1)*(z + 1)*(x + y + z - 2), (x + 1)*(y + 1)*(z + 1) +
                              (x + 1)*(z + 1)*(x + y + z - 2), (x + 1)*(y + 1)*(z + 1) + (x + 1)*(y + 1)*(x + y + z - 2)],
                             [-(1 - x)*(y + 1)*(z + 1) - (y + 1)*(z + 1)*(-x + y + z - 2), (1 - x)*(y + 1)*(z + 1) + (1 - x)
                              * (z + 1)*(-x + y + z - 2), (1 - x)*(y + 1)*(z + 1) + (1 - x)*(y + 1)*(-x + y + z - 2)],
                             [-4*x*(1 - y)*(1 - z), (2 - 2*x**2) *
                              (z - 1), (2 - 2*x**2)*(y - 1)],
                             [2*(1 - y**2)*(1 - z), -2*y*(1 - z) *
                              (2*x + 2), (2*x + 2)*(y**2 - 1)],
                             [-4*x*(1 - z)*(y + 1), (1 - z) *
                              (2 - 2*x**2), (2 - 2*x**2)*(-y - 1)],
                             [-2*(1 - y**2)*(1 - z), -2*y*(1 - z) *
                              (2 - 2*x), (2 - 2*x)*(y**2 - 1)],
                             [-2*(1 - y)*(1 - z**2), (2 - 2*x) *
                              (z**2 - 1), -2*z*(1 - y)*(2 - 2*x)],
                             [2*(1 - y)*(1 - z**2), (2*x + 2) *
                                 (z**2 - 1), -2*z*(1 - y)*(2*x + 2)],
                             [2*(1 - z**2)*(y + 1), (1 - z**2) *
                                 (2*x + 2), -2*z*(2*x + 2)*(y + 1)],
                             [-2*(1 - z**2)*(y + 1), (1 - z**2) *
                              (2 - 2*x), -2*z*(2 - 2*x)*(y + 1)],
                             [-4*x*(1 - y)*(z + 1), (2 - 2*x**2) *
                              (-z - 1), (1 - y)*(2 - 2*x**2)],
                             [2*(1 - y**2)*(z + 1), -2*y*(2*x + 2)
                              * (z + 1), (1 - y**2)*(2*x + 2)],
                             [-4*x*(y + 1)*(z + 1), (2 - 2*x**2) *
                              (z + 1), (2 - 2*x**2)*(y + 1)],
                             [-2*(1 - y**2)*(z + 1), -2*y*(2 - 2*x)*(z + 1), (1 - y**2)*(2 - 2*x)]])
