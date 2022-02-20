"""
*********************************************
Creation of 2D elements
*********************************************

Coords and gdls (degrees of freedom) are given by a Numpy ndarray matrix.

In the coordinate matrix, each row represents a node of the element and each column a dimension.
For example, a 2D triangular element of 3 nodes must have a 3x2 matrix.


In the gdls matrix, each row represents a variable of the node and each column a node.
For example, a 2D triangular element of 3 nodes and 2 variables per node (plane stress) must have a 2x3 gdls matrix.


In this example several characteristics are tested:

- The element creation and transformations are tested over 2 triangular (2D), 2 quadrilateral (2D) and 2 lineal (1D) elements.
- The isBetwwen methods gives the oportunity to check if a given set of points is are inside a element.
- The inverseMapping method allows to give a set of global coordinates and convert them to natural coordinates.
- The jacobian graph allows to verigfy the numerical estability of the element

Coordinate trasformation and shape functions
############################################

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example1_1.png
   :align: center

   Serendipity (8 nodes 2D) element coordinate transformation and shape functions.

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example1_2.png
   :align: center

   Quadrilateral (4 nodes 2D) element coordinate transformation and shape functions.

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example1_3.png
   :align: center

   Triangular (6 nodes 2D) element coordinate transformation and shape functions.

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example1_4.png
   :align: center

   Triangular (3 nodes 2D) element coordinate transformation and shape functions.

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example1_5.png
   :align: center

   Line (3 nodes 1D) element coordinate transformation and shape functions.

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example1_6.png
   :align: center

   Line (2 nodes 2D) element coordinate transformation and shape functions.

Point inside element (works in all dimensions)
##############################################

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example1_7.png
   :align: center

   Test if given points are inside the element.

Jacobian graphs
###############

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example1_8.png
   :align: center

   Serendipity (8 nodes 2D) element jacobian graph.

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example1_9.png
   :align: center

   Quadrilateral (4 nodes 2D) element jacobian graph.

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example1_10.png
   :align: center

   Triangular (6 nodes 2D) element jacobian graph.

.. figure:: https://raw.githubusercontent.com/ZibraMax/FEM/master/Examples/examples_results/example1_11.png
   :align: center

   Triangular (3 nodes 2D) element jacobian graph.

Code
####

.. code:: python

    from FEM.Elements.E2D.Serendipity import Serendipity
    from FEM.Elements.E2D.Quadrilateral import Quadrilateral
    from FEM.Elements.E2D.LTriangular import LTriangular
    from FEM.Elements.E2D.QTriangular import QTriangular

    from FEM.Elements.E1D.LinealElement import LinealElement
    from FEM.Elements.E1D.QuadraticElement import QuadraticElement

    elements = []

    coords = [[1, 1], [3, 2], [3.5, 3], [0.5, 4], [2, 1.5],
              [3.25, 2.5], [(3.5+.5)/2, 3.5], [(0.5+1)/2, 5/2]]
    gdl = np.array([[1, 2, 3, 4, 5, 6, 7, 8]])
    eRS = Serendipity(coords, gdl)
    eRS2 = Serendipity(coords, gdl)
    elements.append(eRS)

    coords = [[1, 1], [3, 2], [3.5, 3], [0.5, 4]]
    gdl = np.array([[1, 2, 3, 4]])
    eRC = Quadrilateral(coords, gdl)
    elements.append(eRC)

    # [[-1.5, -1], [1, 0], [0, 5], [-0.25, -0.5], [0.5, 2.5], [-0.75, 2.0]]
    coords = [[-1.5, -1], [1, 0], [0, 5]]
    for i in range(len(coords)-1):
        coords += [[(coords[i][0]+coords[i+1][0])/2,
                    (coords[i][1]+coords[i+1][1])/2]]
    coords += [[(coords[i+1][0]+coords[0][0])/2,
                (coords[i+1][1]+coords[0][1])/2]]
    gdl = np.array([[1, 2, 3, 4, 5, 6]])
    eTC = QTriangular(coords, gdl)
    elements.append(eTC)

    coords = [[1.4, 1.5], [3.1, 2.4], [2, 3]]
    gdl = np.array([[1, 2, 3]])
    eTL = LTriangular(coords, gdl)
    elements.append(eTL)

    coords = [[3], [4], [5]]
    gdl = np.array([[3, 4, 5]])
    e = QuadraticElement(coords, gdl, 3)
    elements.append(e)

    coords = [[3], [5]]
    gdl = np.array([[3, 5]])
    e = LinealElement(coords, gdl, 3)
    elements.append(e)

    # Coordinate transformation
    for i, e in enumerate(elements):
        e.draw()
        plt.show()

    # Point inside element
    p_test = np.array([1.25, 5.32, 3.1, 3.5, 5.0])
    result = e.isInside(p_test)
    e.draw()
    plt.gca().plot(p_test[result], [0.0] *
                   len(p_test[result]), 'o', c='g', label='Inside')
    plt.gca().plot(p_test[np.invert(result)], [0.0] *
                   len(p_test[np.invert(result)]), 'o', c='r', label='Not Inside')
    plt.legend()
    plt.legend()
    plt.show()

    # Inverse Mapping
    z = eTL.inverseMapping(np.array([[1.3, 2.5, 3.5], [1.5, 2.6, 8.5]]))
    # z = eTL.inverseMapping(np.array([[1.3],[1.5]]))
    print(z)
    # print(eTL.isInside(np.array([3.5,2.5])))

    # Jacobian Graphs
    for i, e in enumerate(elements):
        e.jacobianGraph()
        plt.show()

"""

if __name__ == '__main__':

    # Creación de elementos
    import numpy as np
    import matplotlib.pyplot as plt
    try:
        from FEM.Elements.E2D.Serendipity import Serendipity
        from FEM.Elements.E2D.Quadrilateral import Quadrilateral
        from FEM.Elements.E2D.LTriangular import LTriangular
        from FEM.Elements.E2D.QTriangular import QTriangular

        from FEM.Elements.E1D.LinealElement import LinealElement
        from FEM.Elements.E1D.QuadraticElement import QuadraticElement
    except Exception as e:
        import os
        import inspect
        import sys
        currentdir = os.path.dirname(os.path.abspath(
            inspect.getfile(inspect.currentframe())))
        parentdir = os.path.dirname(currentdir)
        sys.path.insert(0, parentdir)
        from FEM.Elements.E2D.Serendipity import Serendipity
        from FEM.Elements.E2D.Quadrilateral import Quadrilateral
        from FEM.Elements.E2D.LTriangular import LTriangular
        from FEM.Elements.E2D.QTriangular import QTriangular

        from FEM.Elements.E1D.LinealElement import LinealElement
        from FEM.Elements.E1D.QuadraticElement import QuadraticElement

    elements = []

    coords = [[1, 1], [3, 2], [3.5, 3], [0.5, 4], [2, 1.5],
              [3.25, 2.5], [(3.5+.5)/2, 3.5], [(0.5+1)/2, 5/2]]
    gdl = np.array([[1, 2, 3, 4, 5, 6, 7, 8]])
    eRS = Serendipity(coords, gdl)
    eRS2 = Serendipity(coords, gdl)
    elements.append(eRS)

    coords = [[1, 1], [3, 2], [3.5, 3], [0.5, 4]]
    gdl = np.array([[1, 2, 3, 4]])
    eRC = Quadrilateral(coords, gdl)
    elements.append(eRC)

    # [[-1.5, -1], [1, 0], [0, 5], [-0.25, -0.5], [0.5, 2.5], [-0.75, 2.0]]
    coords = [[-1.5, -1], [1, 0], [0, 5]]
    for i in range(len(coords)-1):
        coords += [[(coords[i][0]+coords[i+1][0])/2,
                    (coords[i][1]+coords[i+1][1])/2]]
    coords += [[(coords[i+1][0]+coords[0][0])/2,
                (coords[i+1][1]+coords[0][1])/2]]
    gdl = np.array([[1, 2, 3, 4, 5, 6]])
    eTC = QTriangular(coords, gdl)
    elements.append(eTC)

    coords = [[1.4, 1.5], [3.1, 2.4], [2, 3]]
    gdl = np.array([[1, 2, 3]])
    eTL = LTriangular(coords, gdl)
    elements.append(eTL)

    coords = [[3], [4], [5]]
    gdl = np.array([[3, 4, 5]])
    e = QuadraticElement(coords, gdl, 3)
    elements.append(e)

    coords = [[3], [5]]
    gdl = np.array([[3, 5]])
    e = LinealElement(coords, gdl, 3)
    elements.append(e)

    # Coordinate transformation
    for i, e in enumerate(elements):
        e.draw()
        plt.savefig(f'Examples/examples_results/example1_{i+1}.png')
        plt.close()

    # Point inside element
    p_test = np.array([1.25, 5.32, 3.1, 3.5, 5.0])
    result = e.isInside(p_test)
    e.draw()
    plt.gca().plot(p_test[result], [0.0] *
                   len(p_test[result]), 'o', c='g', label='Inside')
    plt.gca().plot(p_test[np.invert(result)], [0.0] *
                   len(p_test[np.invert(result)]), 'o', c='r', label='Not Inside')
    plt.legend()
    plt.legend()
    plt.savefig(f'Examples/examples_results/example1_{7}.png')
    plt.close()

    # Inverse Mapping
    z = eTL.inverseMapping(np.array([[1.3, 2.5, 3.5], [1.5, 2.6, 8.5]]))
    # z = eTL.inverseMapping(np.array([[1.3],[1.5]]))
    print(z)
    # print(eTL.isInside(np.array([3.5,2.5])))

    # Jacobian Graphs
    for i, e in enumerate(elements):
        e.jacobianGraph()
        plt.savefig(f'Examples/examples_results/example1_{8+i}.png')
        plt.close()
