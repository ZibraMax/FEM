"""
**********************************************
Numerical validation of the Elasticity3D Class
**********************************************

Tests
#####

1 tests where made for validating the Elasticity2D class.

All of these 1 test are based in the same procedure. The objective is to compare the solution of the displacements profile in a cantilever beam.


Geometry
########

The input geometry has 10000 elements.

The beam geometrical properties are:

- b = 0.3 m
- h = 0.6 m
- L = 3.5 m
- E = 21000 KPa
- :math:`\\gamma=23.54\\frac{kN}{m^3}`
- :math:`v=0.2`

All tests are compared whit the analitycal solution of the beam deflections.

The results are obtained using the profile option of the package. The location of all profiles are in the centroid of the beam alog the y axis.

The Plane Stress formulation include shear effects. So, It is needed to obtain analytical solutions that includes such effects. The easyest way is to obtain the analytical solution using the virtual work statement.

Virtual work procedure for a cantilever beam
############################################

The virtual work statement is useful for calculate the dispacements in a given point. The displacements with felxion and shear effects can be calculated as:

.. math::

    \\hat{p}U=\\int_{0}^{L} \\frac{\\hat{M}M}{EI} dx + \\int_{0}^{L} \\frac{\\hat{V}V}{kAG} dx

where:

- :math:`\\hat{p}` is point load applied in the place where the displacement want to be calculated. The load magnitude is 1.
- :math:`\\hat{M}` is the moment diagram produced by the load :math:`\\hat{p}`.
- :math:`\\hat{V}` is the shear diagram produced by the load :math:`\\hat{p}`.

For calculate the displacements in any point of the element lenght, the load :math:`\\hat{p}` have to be applied in a distance :math:`z`.

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/resources/VW.png
   :align: center
   :alt: Virtual work statement

   Moment and shear diagram for virtual load.

Since the diagrams are only different from 0 in length z, the virtual work statement can be rewritten as:

.. math::

    \\hat{p}U=\\int_{0}^{z} \\frac{\\hat{M}M}{EI} dx + \\int_{0}^{z} \\frac{\\hat{V}V}{kAG} dx

:math:`M` and :math:`V` are case dependent an these expresions will be developed for each test case.

"""
from FEM.Geometry.Geometry import Geometry3D
from FEM.Elasticity3D import Elasticity
import matplotlib.pyplot as plt
import unittest
import numpy as np

TOL = 0.1


class TestElasticity3D(unittest.TestCase):
    """Plane Stress tests
    1. Uniform load applied as an uniform load in all elements.
    """

    def test_cantilever_beam_uniform(self):
        """Test if the FEM solution matches the analytical solution:

        For the load condition, the moment an shear diagram can be described by the following equations:

        .. math::

            M=WLx-\\frac{WL^2}{2}-\\frac{Wx^2}{2}

        .. math::

            V=W(L-x)

        Solving the virtual work statement, the analitycal displacement can be calculated as:

        .. math::

            U = \\frac{Wx^2}{24EI}\\left(x^2+6L^2-4Lx\\right)+\\frac{W}{kAG}\\left(Lx-\\frac{x^2}{2}\\right)

        .. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/resources/results/test_cantilever_beam_uniform_1.png
            :alt: Virtual work statement

            Test result: FEM solution compared whit analytical solution

        """
        E = 21000000.0
        v = 0.2
        h = 0.6
        b = 0.3
        L = 3.5
        gamma = 23.54

        _a = L
        _b = h
        _c = b

        nx = 100
        ny = 10
        nz = 10

        dx = _a/nx
        dy = _b/ny
        dz = _c/nz

        coords = []

        for i in range(nx+1):
            x = i*dx
            for j in range(ny+1):
                y = j*dy
                for k in range(nz+1):
                    z = k*dz
                    coords += [[x, y, z]]

        dicc = []

        def node(i, j, k): return i*(ny+1)*(nz+1)+j*(nz+1)+k

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    node1 = node(i, j, k)
                    node2 = node(i+1, j, k)
                    node3 = node(i+1, j+1, k)
                    node4 = node(i, j+1, k)

                    node5 = node(i, j, k+1)
                    node6 = node(i+1, j, k+1)
                    node7 = node(i+1, j+1, k+1)
                    node8 = node(i, j+1, k+1)

                    dicc += [[node1, node2, node3, node4,
                              node5, node6, node7, node8]]

        def fy(x): return -gamma

        geometria = Geometry3D(
            dicc, coords, ["B1V"]*len(dicc), nvn=3, fast=True)
        cbe = []
        for i in range(len(coords)):
            if 0.0 == coords[i][0]:
                cbe += [[i*3, 0.0]]
                cbe += [[i*3+1, 0.0]]
                cbe += [[i*3+2, 0.0]]
        geometria.cbe = cbe

        O = Elasticity(geometria, E, v, gamma, fy=fy,
                       verbose=True, name='3D test_cantilever_beam_uniform')
        O.solve()

        A = b*h
        kk = 5/6
        G = E/(2*(1+v))
        W = gamma*A
        I = b*h**3/12
        print(len(O.geometry.detectBorderElements()), len(O.elements))

        O.exportJSON("3D_BEAM.json")

        def g(x): return W*x**2/24/E/I * \
            (x**2+6*L**2-4*L*x)+W*(L*x-x**2/2)/kk/A/G
        plt.close()
        anal = g(L)
        err = np.max(np.abs(O.U)) - anal
        errores = np.max(np.abs(err/anal))

        # plt.plot(X, U, '--*', color='black', label='FEM')
        # plt.plot(X, anal, '-^', color='black', label='Analytical')
        # plt.grid()
        # plt.legend()
        # plt.xlabel('x [m]')
        # plt.ylabel(r'$\Delta [m]$')
        # plt.savefig('Test/resources/results/test_cantilever_beam_uniform_1.png')
        # plt.close()
        self.assertTrue(errores < TOL)


if __name__ == '__main__':
    unittest.main()
