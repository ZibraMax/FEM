"""
********************************************
Geometry creation using triangles. Torsion2D
********************************************

An square shape section when a unitary rotation is applied.

The input mesh is created using the Delaunay class using second order elements.

Geomery
#######

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example3_geometry.png
   :align: center

   Input geometry created using Delaunay triangulations.

Result
######

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example3.png
   :align: center

   Analysis results.

Code
####

.. code:: python

    import matplotlib.pyplot as plt
    from FEM.Torsion2D import Torsion2D
    from FEM.Geometry import Geometry2D, Delaunay

    a = 1
    b = 1
    E = 200000
    v = 0.27
    G = E/(2*(1+v))*0+1
    phi = 1

    # Use these lines to create the input geometry and file
    # coords = np.array([[0, 0], [1, 0], [1, 1], [0, 1.0]])
    # params = Delaunay._strdelaunay(a=0.0001, o=2)
    # geo = Delaunay(coords, params)
    # geo.exportJSON('Examples/Mesh_tests/Square_torsion.json')

    geometria = Geometry2D.importJSON(
        'Examples/Mesh_tests/Square_torsion.json')
    geometria.show()
    plt.show()
    O = Torsion2D(geometria, G, phi, verbose=True)
    O.solve()
    plt.show()

"""

if __name__ == '__main__':
    try:
        import matplotlib.pyplot as plt
        from FEM.Torsion2D import Torsion2D
        from FEM.Geometry import Geometry2D, Delaunay
    except Exception as e:
        import os
        import inspect
        import sys
        currentdir = os.path.dirname(os.path.abspath(
            inspect.getfile(inspect.currentframe())))
        parentdir = os.path.dirname(currentdir)
        sys.path.insert(0, parentdir)
        import matplotlib.pyplot as plt
        from FEM.Torsion2D import Torsion2D
        from FEM.Geometry import Geometry2D, Delaunay

    a = 1
    b = 1
    E = 200000
    v = 0.27
    G = E/(2*(1+v))*0+1
    phi = 1
    # import numpy as np
    # coords = np.array([[0, 0], [1, 0], [1, 1], [0, 1.0]])
    # params = Delaunay._strdelaunay(a=0.001, o=2)
    # geo = Delaunay(coords, params)
    # geo.exportJSON('Examples/Mesh_tests/Square_torsion.json')
    geometria = Geometry2D.importJSON(
        'Examples/Mesh_tests/Square_torsion.json')
    geometria.show()
    plt.savefig('Examples/examples_results/example3_geometry.png')
    plt.show()
    O = Torsion2D(geometria, G, phi, verbose=True)
    O.solve()
    plt.savefig('Examples/examples_results/example3.png')
    plt.show()
