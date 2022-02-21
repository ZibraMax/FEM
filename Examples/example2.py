"""
********************************************
Geometry creation using triangles. Torsion2D
********************************************

An I shape section when a unitary rotation is applied.

The input mesh is created using the Delaunay class using second order elements.

Geomery
#######

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example2_geometry.png
   :align: center

   Input geometry created using Delaunay triangulations.

Result
######

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example2.png
   :align: center

   Analysis results.

Code
####

.. code:: python

    import matplotlib.pyplot as plt
    from FEM.Torsion2D import Torsion2D
    from FEM.Geometry import Delaunay, Geometry2D

    a = 0.3
    b = 0.3
    tw = 0.05
    tf = 0.05
    E = 200000
    v = 0.27
    G = E / (2 * (1 + v))
    phi = 1
    # Geometry perimeter
    vertices = [
        [0, 0],
        [a, 0],
        [a, tf],
        [a / 2 + tw / 2, tf],
        [a / 2 + tw / 2, tf + b],
        [a, tf + b],
        [a, 2 * tf + b],
        [0, 2 * tf + b],
        [0, tf + b],
        [a / 2 - tw / 2, tf + b],
        [a / 2 - tw / 2, tf],
        [0, tf],
    ]
    # Input string to the delaunay triangulation.
    # o=2 second order elements.
    # a=0.00003 Maximum element area.
    params = Delaunay._strdelaunay(constrained=True, delaunay=True,
                                   a=0.00003, o=2)
    # Mesh reation.
    geometria = Delaunay(vertices, params)

    # geometria.exportJSON('Examples/Mesh_tests/I_test.json') # Save the mesh to a file.
    # geometria = Geometry2D.importJSON('Examples/Mesh_tests/I_test.json') # Load mesh from file.
    
    # Show the geometry.
    geometria.show()
    plt.show()

    # Create the Torsion2D analysis.
    O = Torsion2D(geometria, G, phi)
    O.solve() # All finite element steps.
    plt.show()

"""

if __name__ == '__main__':

    # Creación de elementos
    try:
        import matplotlib.pyplot as plt
        from FEM.Torsion2D import Torsion2D
        from FEM.Geometry import Delaunay, Geometry2D
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
        from FEM.Geometry import Delaunay, Geometry2D

    a = 0.3
    b = 0.3
    tw = 0.05
    tf = 0.05
    E = 200000
    v = 0.27
    G = E / (2 * (1 + v))
    phi = 1
    vertices = [
        [0, 0],
        [a, 0],
        [a, tf],
        [a / 2 + tw / 2, tf],
        [a / 2 + tw / 2, tf + b],
        [a, tf + b],
        [a, 2 * tf + b],
        [0, 2 * tf + b],
        [0, tf + b],
        [a / 2 - tw / 2, tf + b],
        [a / 2 - tw / 2, tf],
        [0, tf],
    ]
    params = Delaunay._strdelaunay(constrained=True, delaunay=True,
                                   a='0.00009', o=2)
    geometria = Delaunay(vertices, params)

    # geometria.exportJSON('Examples/Mesh_tests/I_test.json')
    # geometria = Geometry2D.importJSON('Examples/Mesh_tests/I_test.json')
    geometria.show()
    plt.savefig(f'Examples/examples_results/example2_geometry.png')
    plt.show()
    print(len(geometria.elements))
    O = Torsion2D(geometria, G, phi)
    O.solve()
    plt.savefig(f'Examples/examples_results/example2.png')
    plt.show()
