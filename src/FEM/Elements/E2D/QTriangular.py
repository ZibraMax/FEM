"""
*********************************************
TRIANGULAR SECOND ORDER ELEMENT
*********************************************

Defines the lagrange second order triangular element

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/docs/source/QTriangular.png
   :align: center

   6 nodes triangular element. (Reddy, 2005)


Shape Functions
---------------

.. math::
    \\Psi_i=\\left[\\begin{matrix}\\left(x + y - 0.5\\right) \\left(2.0 x + 2.0 y - 2.0\\right)\\\\2.0 x \\left(x - 0.5\\right)\\\\2.0 y \\left(y - 0.5\\right)\\\\x \\left(- 4.0 x - 4.0 y + 4.0\\right)\\\\4.0 x y\\\\- 4.0 y \\left(x + y - 1.0\\right)\\end{matrix}\\right]

Shape Functions Derivatives
---------------------------

.. math::
    \\frac{\\partial \\Psi_i}{\\partial x_j}=\\left[\\begin{matrix}4.0 x + 4.0 y - 3.0 & 4.0 x + 4.0 y - 3.0\\\\4.0 x - 1.0 & 0\\\\0 & 4.0 y - 1.0\\\\- 8.0 x - 4.0 y + 4.0 & - 4.0 x\\\\4.0 y & 4.0 x\\\\- 4.0 y & - 4.0 x - 8.0 y + 4.0\\end{matrix}\\right]


"""


from .Element2D import Element2D, np
from .TriangularScheme import TriangularScheme
from ..E1D.QuadraticElement import QuadraticElement


class QTriangular(Element2D, TriangularScheme):
    """Creates a lagrangian element of order 2

    Args:
        coords (np.ndarray): Element coordinates matrix
        gdl (np.ndarray): Element gdl matrix
        n (int, optional): Number of Gauss Points. Defaults to 2.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 3, **kargs) -> None:
        """Creates a lagrangian element of order 2

        Args:
            coords (np.ndarray): Element coordinates matrix
            gdl (np.ndarray): Element gdl matrix
            n (int, optional): Number of Gauss Points. Defaults to 3.
        """

        coords = np.array(coords)
        delta = coords[0]-coords[1]
        he1 = np.linalg.norm(delta)
        e1 = QuadraticElement(np.array([[0], [he1*0.5], [he1]]),
                              np.array([[-1, -1, -1]]), border=True)

        delta = coords[2]-coords[1]
        he2 = np.linalg.norm(delta)
        e2 = QuadraticElement(np.array([[0], [he2*0.5], [he2]]),
                              np.array([[-1, -1, -1]]), border=True)

        delta = coords[0]-coords[2]
        he3 = np.linalg.norm(coords[0]-coords[2])
        e3 = QuadraticElement(np.array([[0], [he3*0.5], [he3]]),
                              np.array([[-1, -1, -1]]), border=True)
        self.borders = [e1, e2, e3]

        _coords = np.array([coords[i] for i in range(3)])
        TriangularScheme.__init__(self, n, **kargs)
        Element2D.__init__(self, coords, _coords, gdl, **kargs)

    def psis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """
        return np.array([
            2.0*(z[0]+z[1]-1.0)*(z[0]+z[1]-0.5),
            2.0*z[0]*(z[0]-0.5),
            2.0*z[1]*(z[1]-0.5),
            -4.0*(z[0]+z[1]-1.0)*(z[0]),
            4.0*z[0]*z[1],
            -4.0*z[1]*(z[0]+z[1]-1.0)]).T

    def dpsis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """

        return np.array([
            [4.0*z[0]+4.0*z[1]-3.0, 4.0*z[1]+4.0*z[0]-3.0],
            [4.0*z[0]-1.0, 0*z[0]],
            [0*z[0], 4.0*z[1]-1.0],
            [-8.0*z[0]-4.0*(z[1]-1.0), -4.0*z[0]],
            [4.0*z[1], 4.0*z[0]],
            [-4.0*z[1], -8.0*z[1]-4.0*z[0]+4.0]
        ])
