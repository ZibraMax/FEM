"""
********************************
Ordinary diferential equation 1D
********************************

The current file ust the general formulation of the EDO1D Class using custom function coeficients.

Geomery
#######

.. math::

    a(x)\\frac{d^2u}{dx^2}+c(x)u=f(x)
.. math::

    a(x)=x^2-2
.. math::

    c(x)=x-3
.. math::

    f(x)=6(x^2-2)+(x-3)(3x^2)

.. math::

    U|_0=0
.. math::

    U|_L=3
.. math::

    L=1

Result
######

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example4.png
   :align: center

   Analysis results.

Code
####

.. code:: python

    import matplotlib.pyplot as plt
    from FEM.EDO1D import EDO1D
    from FEM.Geometry import Lineal

    # Define functions of the diferential equation
    def a(x): return (x[0]**2-2)
    def c(x): return (x[0]-3)
    def f(x): return (x[0]**2-2)*6+(x[0]-3)*(3*x[0]**2)

    # Define border conditions. List of border conditions. In the first node, value=0.0, in the last node, value = 3.0
    cbe = [[0, 0.0], [-1, 3.0]]

    lenght = 1
    n = 500
    o = 2
    geometria = Lineal(lenght, n, o)
    O = EDO1D(geometria, a, c, f)
    O.cbe = cbe
    O.solve()
    plt.show()
"""

if __name__ == '__main__':

    import matplotlib.pyplot as plt
    from FEM.Elasticity2D import PlaneStressSparse
    from FEM.Geometry import Geometry2D
    from FEM.Utils import enmalladoFernando

    a = 5
    u0 = 0.001

    E = 2.1*10**6
    v = 0.2
    t = 0.5

    coords, dicc = enmalladoFernando(a, a, 30, 30)
    geo = Geometry2D(dicc, coords, ['C2V']*len(dicc))
    geo.generateRegionFromCoords([0.0, 0.0], [a, 0.0])
    geo.generateRegionFromCoords([a, 0.0], [a, a])
    geo.generateRegionFromCoords([a, a], [0.0, a])
    geo.generateRegionFromCoords([0.0, a], [0.0, 0.0])
    geo.exportJSON("Examples/Mesh_tests/rect.json")

    geo = Geometry2D(dicc, coords, ['C2V']*len(dicc), nvn=2)
    geo.generateRegionFromCoords([0.0, 0.0], [a, 0.0])
    geo.generateRegionFromCoords([a, 0.0], [a, a])
    geo.generateRegionFromCoords([a, a], [0.0, a])
    geo.generateRegionFromCoords([0.0, a], [0.0, 0.0])
    geo.exportJSON("Examples/Mesh_tests/rect2.json")

    geometria = Geometry2D.importJSON(
        "Examples/Mesh_tests/rect.json", fast=True)
    O = PlaneStressSparse(geometria, E, v, t)
    cbe = O.geometry.cbFromRegion(3, 0, 1)
    cbe += O.geometry.cbFromRegion(3, 0, 2)
    cbe += O.geometry.cbFromRegion(1, u0, 1)
    O.cbe = cbe
    O.solve()
    plt.show()
