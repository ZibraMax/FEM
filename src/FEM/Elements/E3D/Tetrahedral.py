"""
*********************************************
TETRAHEDRAL ELEMENTS
*********************************************

Defines the lagrangian first order and second order tetrahedral elements

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/docs/source/Tetrahedrals.png
   :align: center

   Tetrahedral (4 nodes) and TetrahedralO2 (10 nodes). (Reddy, 2005)

Tetrahedral
###########

First order 4 node tetrahedral element

Shape Functions
---------------

.. math::
    \\Psi_i=\\left[\\begin{matrix}- x - y - z + 1\\\\x\\\\y\\\\z\\end{matrix}\\right]

Shape Functions Derivatives
---------------------------

.. math::
    \\frac{\\partial \\Psi_i}{\\partial x_j}=\\left[\\begin{matrix}-1 & -1 & -1\\\\1 & 0 & 0\\\\0 & 1 & 0\\\\0 & 0 & 1\\end{matrix}\\right]

TetrahedralO2
#############

Second order 10 node tetrahedral element

Shape Functions
---------------

.. math::
    \\Psi_i=\\left[\\begin{matrix}\\left(- 2 x - 2 y - 2 z + 1\\right) \\left(- x - y - z + 1\\right)\\\\x \\left(2 x - 1\\right)\\\\y \\left(2 y - 1\\right)\\\\z \\left(2 z - 1\\right)\\\\x \\left(- 4 x - 4 y - 4 z + 4\\right)\\\\4 x y\\\\4 y \\left(- x - y - z + 1\\right)\\\\z \\left(- 4 x - 4 y - 4 z + 4\\right)\\\\4 x z\\\\4 y z\\end{matrix}\\right]

Shape Functions Derivatives
---------------------------

.. math::
    \\frac{\\partial \\Psi_i}{\\partial x_j}=\\left[\\begin{matrix}4 x + 4 y + 4 z - 3 & 4 x + 4 y + 4 z - 3 & 4 x + 4 y + 4 z - 3\\\\4 x - 1 & 0 & 0\\\\0 & 4 y - 1 & 0\\\\0 & 0 & 4 z - 1\\\\- 8 x - 4 y - 4 z + 4 & - 4 x & - 4 x\\\\4 y & 4 x & 0\\\\- 4 y & - 4 x - 8 y - 4 z + 4 & - 4 y\\\\- 4 z & - 4 z & - 4 x - 4 y - 8 z + 4\\\\4 z & 0 & 4 x\\\\0 & 4 z & 4 y\\end{matrix}\\right]

"""

from ..E2D.LTriangular import LTriangular
from ..E2D.QTriangular import QTriangular
from .Element3D import Element3D, np
from .TetrahedralScheme import TetrahedralScheme


class Tetrahedral(Element3D, TetrahedralScheme):
    """Creates a 3D tetrahedral element

    Args:
        coords (np.ndarray): Node coordinates matrix
        gdl (np.ndarray): Degrees of freedom matrix
        n (int, optional): Number of gauss points used for integration. Defaults to 3.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 3, **kargs) -> None:
        """Creates a 3D tetrahedral element

        Args:
            coords (np.ndarray): Node coordinates matrix
            gdl (np.ndarray): Degrees of freedom matrix
            n (int, optional): Number of gauss points used for integration. Defaults to 3.
        """

        coords = np.array(coords)
        self.faces = [
            [0, 1, 3],
            [1, 2, 3],
            [0, 3, 2],
            [0, 2, 1]]
        self.face_element = LTriangular

        TetrahedralScheme.__init__(self, n, **kargs)
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
        L1 = 1-x-y-z
        L2 = x
        L3 = y
        L4 = z
        return np.array(
            [L1, L2, L3, L4]).T

    def dpsis(self, _z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """
        x = _z[0]
        kernell = x-x

        return np.array([
            [-1.0+kernell, -1.0+kernell, -1.0+kernell],
            [1.0+kernell, 0.0+kernell, 0.0+kernell],
            [0.0+kernell, 1.0+kernell, 0.0+kernell],
            [0.0+kernell, 0.0+kernell, 1.0+kernell]])


class TetrahedralO2(Element3D, TetrahedralScheme):
    """Creates a 3D second order tetrahedral element

    Args:
        coords (np.ndarray): Node coordinates matrix
        gdl (np.ndarray): Degrees of freedom matrix
        n (int, optional): Number of gauss points used for integration. Defaults to 3.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 3, **kargs) -> None:
        """Creates a 3D second order tetrahedral element

        Args:
            coords (np.ndarray): Node coordinates matrix
            gdl (np.ndarray): Degrees of freedom matrix
            n (int, optional): Number of gauss points used for integration. Defaults to 3.
        """

        coords = np.array(coords)
        self.faces = [
            [0, 1, 3, 4, 8, 7],
            [1, 2, 3, 5, 9, 8],
            [0, 3, 2, 7, 9, 6],
            [0, 2, 1, 6, 5, 4]]
        self.face_element = QTriangular

        TetrahedralScheme.__init__(self, n, **kargs)
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
        L1 = 1-x-y-z
        L2 = x
        L3 = y
        L4 = z
        return np.array([
            L1*(2*L1-1),
            L2*(2*L2-1),
            L3*(2*L3-1),
            L4*(2*L4-1),
            4*L1*L2,
            4*L2*L3,
            4*L3*L1,
            4*L1*L4,
            4*L2*L4,
            4*L3*L4]).T

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

        return np.array([
            [4*x + 4*y + 4*z - 3, 4*x + 4*y + 4*z - 3, 4*x + 4*y + 4*z - 3],
            [4*x - 1, 0, 0],
            [0, 4*y - 1, 0],
            [0, 0, 4*z - 1],
            [-8*x - 4*y - 4*z + 4, -4*x, -4*x],
            [4*y, 4*x, 0],
            [-4*y, -4*x - 8*y - 4*z + 4, -4*y],
            [-4*z, -4*z, -4*x - 4*y - 8*z + 4],
            [4*z, 0, 4*x],
            [0, 4*z, 4*y]])
