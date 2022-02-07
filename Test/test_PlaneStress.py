"""
*********************************************
Numerical validation of the PlaneStress Class
*********************************************

Tests
#####

5 tests where made for validating the PlaneStress class.

All of these 5 test are based in the same procedure. The objective is to compare the solution of the displacements profile in a cantilever beam.


Geometry
########

The input geometry has 400 elements, 10 in y direction, 40 in x directions. The input geometry file can be generated using TestGeometry class.

The beam geometrical properties are:

- b = 0.3 m
- h = 0.5 m
- L = 2 m
- E = 20000 KPa
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
   :scale: 80 %
   :align: center
   :alt: Virtual work statement

   Moment and shear diagram for virtual load.

Since the diagrams are only different from 0 in length z, the virtual work statement can be rewritten as:

.. math::

    \\hat{p}U=\\int_{0}^{z} \\frac{\\hat{M}M}{EI} dx + \\int_{0}^{z} \\frac{\\hat{V}V}{kAG} dx

:math:`M` and :math:`V` are case dependent an these expresions will be developed for each test case.

"""

from FEM.Geometry.Geometry import Geometry
from FEM.Elasticity2D import PlaneStress
import matplotlib.pyplot as plt
import unittest
import os
import numpy as np
import sys
import inspect

FILENAME = 'Test/resources/beam_ws.msh'
TOL = 0.1


class TestPlaneStress(unittest.TestCase):
    """Plane Stress tests
    1. Uniform load applied as an uniform load in all elements.
    2. Uniform load applied as a external load in the superior elements of the beam.
    3. Point load applied in the beam lenght. The load is applied as a shear stress load.
    4. Point load applied in the beam lenght. The load is applied as point load in the corner node.
    5. Triangular load applied in the superior elements of the beam.
    """

    def test_cantilever_beam_uniform_1(self):
        """Test if the FEM solution matches the next analytical solution:

        For the load condition, the moment an shear diagram can be described by the following equations:

        .. math::

            M=WLx-\\frac{WL^2}{2}-\\frac{Wx^2}{2}

        .. math::

            V=W(L-x)

        Solving the virtual work statement, the analitycal displacement can be calculated as:

        .. math::

            U = \\frac{Wx^2}{24EI}\\left(x^2+6L^2-4Lx\\right)+\\frac{W}{kAG}\\left(Lx-\\frac{x^2}{2}\\right)

        .. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/resources/results/test_cantilever_beam_uniform_1.png
            :scale: 50 %
            :alt: Virtual work statement

            Test result: FEM solution compared whit analytical solution

        """
        geometry = Geometry.loadmsh(FILENAME)
        cb = geometry.cbFromSegment(3, 0.0, 1)
        cb += geometry.cbFromSegment(3, 0.0, 2)
        geometry.cbe = cb
        geometry.maskFromSegments()
        self.assertTrue(len(geometry.cbe) > 0)
        h = 0.5
        L = 2.0
        E = 20000
        v = 0.2
        b = 0.3
        I = b*h**3/12
        A = b*h
        k = 5/6
        G = E/(2*(1+v))
        gamma = 23.54
        W = A*gamma
        O = PlaneStress(geometry, E, v, b, fy=lambda x: -gamma)
        O.solve(plot=False)
        X, U = O.profile([0, 0.5*h], [L, 0.5*h], n=10)
        U = np.array(U)[:, 1].flatten()
        def g(x): return W*x**2/24/E/I*(x**2+6*L**2-4*L*x)+W*(L*x-x**2/2)/k/A/G
        plt.close()
        anal = -g(np.array(X))
        err = U - anal
        err = err[1:]
        errores = np.max(np.abs(err/anal[1:]))

        plt.plot(X, U, '--*', color='black', label='FEM')
        plt.plot(X, anal, '-^', color='black', label='Analytical')
        plt.grid()
        plt.legend()
        plt.xlabel('x [m]')
        plt.ylabel(r'$\Delta [m]$')
        plt.savefig('Test/resources/results/test_cantilever_beam_uniform_1.png')
        plt.close()
        self.assertTrue(errores < TOL)

    def test_cantilever_beam_uniform_2(self):
        """Test if the FEM solution matches the next analytical solution:

        For the load condition, the moment an shear diagram can be described by the following equations:

        .. math::

            M=WLx-\\frac{WL^2}{2}-\\frac{Wx^2}{2}

        .. math::

            V=W(L-x)

        Solving the virtual work statement, the analitycal displacement can be calculated as:

        .. math::

            U = \\frac{Wx^2}{24EI}\\left(x^2+6L^2-4Lx\\right)+\\frac{W}{kAG}\\left(Lx-\\frac{x^2}{2}\\right)

        .. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/resources/results/test_cantilever_beam_uniform_2.png
            :scale: 50 %
            :alt: Virtual work statement

            FEM solution compared whit analytical solution
        """

        geometry = Geometry.loadmsh(FILENAME)
        cb = geometry.cbFromSegment(3, 0.0, 1)
        cb += geometry.cbFromSegment(3, 0.0, 2)
        geometry.cbe = cb
        geometry.maskFromSegments()
        self.assertTrue(len(geometry.cbe) > 0)
        h = 0.5
        L = 2.0
        E = 20000
        v = 0.2
        b = 0.3
        I = b*h**3/12
        A = b*h
        gamma = 23.54
        k = 5/6
        G = E/(2*(1+v))
        W = A*gamma
        O = PlaneStress(geometry, E, v, b)
        O.geometry.loadOnSegment(2, fy=lambda x: -W)
        O.solve(plot=False)
        X, U = O.profile([0, 0.5*h], [L, 0.5*h], n=10)
        U = np.array(U)[:, 1].flatten()
        def g(x): return W*x**2/24/E/I*(x**2+6*L**2-4*L*x)+W*(L*x-x**2/2)/k/A/G
        plt.close()
        anal = -g(np.array(X))
        err = U - anal
        err = err[1:].flatten()
        errores = np.max(np.abs(err/anal[1:]))

        plt.plot(X, U, '--*', color='black', label='FEM')
        plt.plot(X, anal, '-^', color='black', label='Analytical')
        plt.grid()
        plt.legend()
        plt.xlabel('x [m]')
        plt.ylabel(r'$\Delta [m]$')
        plt.savefig('Test/resources/results/test_cantilever_beam_uniform_2.png')
        plt.close()
        self.assertTrue(errores < TOL)

    def test_cantilever_beam_point_1(self):
        """Test if the FEM solution matches the next analytical solution:

        For the load condition, the moment an shear diagram can be described by the following equations:

        .. math::

            M=PL\\left(1-\\frac{x}{L}\\right)

        .. math::

            V=P

        Solving the virtual work statement, the analitycal displacement can be calculated as:

        .. math::

            U = \\frac{P^2}{6EI}\\left(3L-x\\right)+\\frac{Px}{kAG}

        .. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/resources/results/test_cantilever_beam_point_1.png
            :scale: 50 %
            :alt: Virtual work statement

            FEM solution compared whit analytical solution

        """
        h = 0.5
        L = 2.0
        E = 20000
        v = 0.2
        b = 0.3
        I = b*h**3/12
        A = b*h
        gamma = 23.54
        W = A*gamma
        G = E/(2*(1+v))
        k = 5/6

        P = 1

        geometry = Geometry.loadmsh(FILENAME)
        cb = geometry.cbFromSegment(3, 0.0, 1)
        cb += geometry.cbFromSegment(3, 0.0, 2)
        geometry.cbe = cb

        geometry.maskFromSegments()
        self.assertTrue(len(geometry.cbe) > 0)

        O = PlaneStress(geometry, E, v, b)
        O.geometry.loadOnSegment(1, fy=lambda x: -P/A*b)
        O.solve(plot=False)
        X, U = O.profile([0, 0.5*h], [L, 0.5*h], n=10)
        U = np.array(U)[:, 1].flatten()
        def g(x): return P*x**2*(3*L-x)/6/E/I+P*x/k/A/G
        plt.close()
        anal = -g(np.array(X))
        err = U - anal
        err = err[1:].flatten()
        errores = np.max(np.abs(err/anal[1:]))

        plt.plot(X, U, '--*', color='black', label='FEM')
        plt.plot(X, anal, '-^', color='black', label='Analytical')
        plt.grid()
        plt.legend()
        plt.xlabel('x [m]')
        plt.ylabel(r'$\Delta [m]$')
        plt.savefig('Test/resources/results/test_cantilever_beam_point_1.png')
        plt.close()
        self.assertTrue(errores < TOL)

    def test_cantilever_beam_point_2(self):
        """Test if the FEM solution matches the next analytical solution:

        For the load condition, the moment an shear diagram can be described by the following equations:

        .. math::

            M=PL\\left(1-\\frac{x}{L}\\right)

        .. math::

            V=P

        Solving the virtual work statement, the analitycal displacement can be calculated as:

        .. math::

            U = \\frac{P^2}{6EI}\\left(3L-x\\right)+\\frac{Px}{kAG}

        .. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/resources/results/test_cantilever_beam_point_2.png
            :scale: 50 %
            :alt: Virtual work statement

            FEM solution compared whit analytical solution

        """
        h = 0.5
        L = 2.0
        E = 20000
        v = 0.2
        b = 0.3
        I = b*h**3/12
        A = b*h
        gamma = 23.54
        W = A*gamma
        G = E/(2*(1+v))
        k = 5/6

        P = 1

        geometry = Geometry.loadmsh(FILENAME)
        cb = geometry.cbFromSegment(3, 0.0, 1)
        cb += geometry.cbFromSegment(3, 0.0, 2)
        cbn = geometry.generateBCFromCoords(L, h, -P, 2)
        geometry.cbe = cb
        geometry.cbn = cbn

        geometry.maskFromSegments()
        self.assertTrue(len(geometry.cbe) > 0)

        O = PlaneStress(geometry, E, v, b)
        O.solve(plot=False)
        X, U = O.profile([0, 0.5*h], [L, 0.5*h], n=10)
        U = np.array(U)[:, 1].flatten()
        def g(x): return P*x**2*(3*L-x)/6/E/I+P*x/k/A/G
        plt.close()
        anal = -g(np.array(X))
        err = U - anal
        err = err[1:].flatten()
        errores = np.max(np.abs(err/anal[1:]))

        plt.plot(X, U, '--*', color='black', label='FEM')
        plt.plot(X, anal, '-^', color='black', label='Analytical')
        plt.grid()
        plt.legend()
        plt.xlabel('x [m]')
        plt.ylabel(r'$\Delta [m]$')
        plt.savefig('Test/resources/results/test_cantilever_beam_point_2.png')
        plt.close()
        self.assertTrue(errores < TOL)

    def test_cantilever_beam_triangular_3(self):
        """Test if the FEM solution matches the next analytical solution:

        For the load condition, the moment an shear diagram can be described by the following equations:

        .. math::

            M=\\frac{-WL^2}{6}+\\frac{LW}{2}x-\\frac{x}{2}\\left(Wx\\left(1-\\frac{x}{L}\\right)\\right)-\\frac{2}{3}x\\left(\\frac{Wx^2}{2L}\\right)

        .. math::

            V=\\frac{Wx^2}{2L}+Wx\\left(1-\\frac{x}{L}\\right)-\\frac{LW}{2}

        Solving the virtual work statement, the analitycal displacement can be calculated as:

        .. math::

            U = \\frac{Wx^2}{120LEI}\\left(10L^3-10L^2x+5Lx^2-x^3\\right)+\\frac{Wx}{6kAGL}\\left(x^2-3Lx+3L^2\\right)

        .. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Test/resources/results/test_cantilever_beam_triangular_3.png
            :scale: 50 %
            :alt: Virtual work statement

            FEM solution compared whit analytical solution

        """
        geometry = Geometry.loadmsh(FILENAME)
        cb = geometry.cbFromSegment(3, 0.0, 1)
        cb += geometry.cbFromSegment(3, 0.0, 2)
        geometry.cbe = cb
        geometry.maskFromSegments()
        self.assertTrue(len(geometry.cbe) > 0)
        h = 0.5
        L = 2.0
        E = 20000
        v = 0.2
        b = 0.3
        I = b*h**3/12
        A = b*h
        gamma = 23.54
        W = A*gamma
        G = E/(2*(1+v))
        k = 5/6
        O = PlaneStress(geometry, E, v, b)
        O.geometry.loadOnSegment(2, fy=lambda x: -W*(x/L))
        O.solve(plot=False)
        X, U = O.profile([0, 0.5*h], [L, 0.5*h], n=10)
        U = np.array(U)[:, 1].flatten()

        def g(x): return W*x**2/120/E/I/L*(10*L**3-10*L**2*x +
                                           5*L*x**2-x**3)+W*x/k/A/G/6/L*(x**2-3*L*x+3*L**2)
        plt.close()
        anal = -g(np.array(X))
        err = U - anal
        err = err[1:].flatten()
        errores = np.max(np.abs(err/anal[1:]))

        plt.plot(X, U, '--*', color='black', label='FEM')
        plt.plot(X, anal, '-^', color='black', label='Analytical')
        plt.grid()
        plt.legend()
        plt.xlabel('x [m]')
        plt.ylabel(r'$\Delta [m]$')
        plt.savefig(
            'Test/resources/results/test_cantilever_beam_triangular_3.png')
        plt.close()
        self.assertTrue(errores < TOL)


if __name__ == '__main__':
    currentdir = os.path.dirname(os.path.abspath(
        inspect.getfile(inspect.currentframe())))
    parentdir = os.path.dirname(currentdir)
    sys.path.insert(0, parentdir)
    from FEM.Geometry.Geometry import Geometry
    from FEM.Elasticity2D import PlaneStress
    unittest.main()
