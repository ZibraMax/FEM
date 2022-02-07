"""
*********************************************
Numerical validation of the Trosion2D Class
*********************************************

Tests
#####

2 tests where made for validating the Torsion2D class. The objective of the test is to calculate the calculating polar moment of inertia (J) of a section.

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/resources/CS.png
   :scale: 80 %
   :align: center
   :alt: Problem statement

Calculating polar moment of inertia (J)
#######################################

Using the formulation of the Torsion2D class, the polar moment of inertia can be calculated as:

.. math::

    J = \\frac{2}{G}\\int_{A}{\\Psi}dA

Where :math:`\\Psi` is the main variable solution of the Torsion2D class.

The integral can be computed by element using the following expresion:"

.. math::

    \\int_{A}{\\Psi}dA = \\sum_{j=1}^{n}I(e_j)

.. math::

    I(e)= \\sum_{i=1}^{N}(U(\\zeta_i)|Jac_i(\\zeta_i)|)W_i

Where:

- :math:`e` is a single area element.
- :math:`|Jac(\\zeta_i)|` is the determinant of the jacobian matrix.
- :math:`U(\\zeta_i)` is the interpolated solution of the main variable.
- :math:`\\zeta_i` is the Gauss-Legendre quadrature point.
- :math:`W_i` is the weight of the i-Gauss point.
- :math:`N` is the number of Gauss points.
- :math:`n` is the number of elements.

Hollow shapes
#############

Due to the boundary conditions of the equation, it is not possible to introduce holes in the geometry.
The holes can be modeled using elements affecting the shear modulus. The hole elements will have relatively close to zero shear modulus, while the material elements will have a higher shear modulus.
In this case, a shear modulus of 80000 is used for material elements and a shear modulus of 1*10^-6 is used for hole elements.

A rotation angle of 1 is selected for all test cases.


"""

from FEM.Geometry.Delaunay import Delaunay
from FEM.Utils.polygonal import giveCoordsCircle
from FEM.Torsion2D import Torsion2D
import matplotlib.path as mpltPath
import unittest
import os
import numpy as np
import sys
import inspect

TOL = 0.01


class TestTorsion2D(unittest.TestCase):
    """Torsion2D Tests
    1. Hollow circle.
    2. Circle.
    """

    def test_inertia_hollow_circle(self):
        """Test if the polar moment of inertia (J) matches the analytical solution for a hollow circle. In this case, the analytical solutions is given by:

        .. math::

            J = \\frac{\\pi}{2}\\left(r_1^4-r_2^4\\right)

        Where:

        - :math:`r_1` is the external radius.
        - :math:`r_2` is the inner radius.
        """

        r = 1.0
        r2 = r/2
        N = 100
        phi = 1.0
        G = 1000.0
        vert, seg = giveCoordsCircle([0, 0], r, n=N)
        vertextra, segextra = giveCoordsCircle([0, 0], r2, n=int(N/2))
        holes = []
        hole = {'center': [-2*r, 2*r],
                'segments': segextra, 'vertices': vertextra}
        holes += [hole]
        params = Delaunay._strdelaunay(
            constrained=True, delaunay=True, a='0.003', o=2)
        geometria = Delaunay(vert, params, nvn=1, holes_dict=holes)
        GG = []
        for centroide in geometria.centroids:
            path = mpltPath.Path(vertextra)
            inside2 = path.contains_points([centroide])
            GG += [G]
            if inside2[0]:
                GG[-1] = 1*10**-6
        geometria.segments = seg
        O = Torsion2D(geometria, GG, phi)
        O.solve()
        integral = 0
        for e in O.elements:
            _, _u = e.giveSolution(domain='gauss-points')
            jac, _ = e.J(e.Z.T)
            detjac = np.linalg.det(jac)
            integral += np.sum(_u*e.W*detjac)
            teorico = np.pi/2*(r**4-r2**4)
        errores = abs((integral*2/G - teorico)/teorico)
        self.assertTrue(errores < TOL)

    def test_inertia_circle(self):
        """Test if the polar moment of inertia (J) matches the analytical solution for circle. In this case, the analytical solutions is given by:

        .. math::

            J = \\frac{\\pi}{2}\\left(r^4\\right)

        Where:

        - :math:`r` is the radius.
        """
        r = 1.0
        N = 100
        phi = 1.0
        G = 1000.0
        vert, _ = giveCoordsCircle([0, 0], r, n=N)
        params = Delaunay._strdelaunay(
            constrained=True, delaunay=True, a='0.003', o=2)
        geometria = Delaunay(vert, params, nvn=1)
        O = Torsion2D(geometria, G, phi)
        O.solve()
        integral = 0
        for e in O.elements:
            _, _u = e.giveSolution(domain='gauss-points')
            jac, _ = e.J(e.Z.T)
            detjac = np.linalg.det(jac)
            integral += np.sum(_u*e.W*detjac)
            teorico = np.pi/2*(r**4)
        errores = abs((integral*2/G - teorico)/teorico)
        self.assertTrue(errores < TOL)


if __name__ == '__main__':
    currentdir = os.path.dirname(os.path.abspath(
        inspect.getfile(inspect.currentframe())))
    parentdir = os.path.dirname(currentdir)
    sys.path.insert(0, parentdir)
    from FEM.Geometry.Delaunay import Delaunay
    from FEM.Torsion2D import Torsion2D
    from FEM.Utils.polygonal import giveCoordsCircle
    unittest.main()
