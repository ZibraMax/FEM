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
    from FEM.EDO1D import EDO1D
    from FEM.Geometry import Lineal

    def a(x): return (x[0]**2-2)
    def c(x): return (x[0]-3)
    def f(x): return (x[0]**2-2)*6+(x[0]-3)*(3*x[0]**2)

    cbe = [[0, 0.0], [-1, 3.0]]

    lenght = 1
    n = 500
    o = 2
    geometria = Lineal(lenght, n, o)
    O = EDO1D(geometria, a, c, f)
    O.cbe = cbe
    O.solve()
    O.exportJSON("Examples/Mesh_tests/Example4.json")
    plt.savefig('Examples/examples_results/example4.png')
    plt.show()
