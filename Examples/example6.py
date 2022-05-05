"""
*********************************************
Test of point inside volumen (3D elements)
*********************************************

This example creates a 3D Brick element and test if a given set of points are inside the element.

Coords and gdls (degrees of freedom) are given by a Numpy ndarray matrix.

In the coordinate matrix, each row represents a node of the element and each column a dimension.
For example, a 2D triangular element of 3 nodes must have a 3x2 matrix.

In the gdls matrix, each row represents a variable of the node and each column a node.
For example, a 2D triangular element of 3 nodes and 2 variables per node (plane stress) must have a 2x3 gdls matrix.


Code
####

.. code:: python

    from FEM.Elements.E3D.Brick import Brick
    import numpy as np

    coords = np.array([[0, 0, 0.0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [
        0, 0, 1.0], [1, 0, 1], [1.5, 1.5, 1.5], [0, 1, 1]])
    coords2 = coords + 1.5

    gdl = np.array([[0, 0, 0, 0, 0, 0, 0, 0]])

    e = Brick(coords=coords, gdl=gdl)
    e2 = Brick(coords=coords2, gdl=gdl)

    domain = e.T(e.domain.T)[0]
    domain2 = e2.T(e.domain.T)[0]

    points = np.array(
        np.array([[-1, 0, 0], [0.5, 0.5, 0.5]]).tolist()+domain.tolist())

    r1 = e.isInside(points)
    r2 = e.isInside(domain2)
    
    print(r1, r2)
    a = 0

"""
if __name__ == '__main__':
    from FEM.Elements.E3D.Brick import Brick
    import numpy as np
    coords = np.array([[0, 0, 0.0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [
        0, 0, 1.0], [1, 0, 1], [1.5, 1.5, 1.5], [0, 1, 1]])
    coords2 = coords + 1.5
    gdl = np.array([[0, 0, 0, 0, 0, 0, 0, 0]])
    e = Brick(coords=coords, gdl=gdl)
    e2 = Brick(coords=coords2, gdl=gdl)
    domain = e.T(e.domain.T)[0]
    domain2 = e2.T(e.domain.T)[0]
    points = np.array(
        np.array([[-1, 0, 0], [0.5, 0.5, 0.5]]).tolist()+domain.tolist())

    r1 = e.isInside(points)
    r2 = e.isInside(domain2)
    print(r1, r2)
    a = 0
