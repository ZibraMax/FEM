if __name__ == '__main__':
    from FEM.Geometry.Region import Region2D
    import scipy
    import json
    import numpy as np
    import matplotlib.pyplot as plt
    from FEM.Elasticity3D import Elasticity
    from FEM.Geometry import Geometry3D, Delaunay
    from FEM.Solvers.Lineal import LinealEigen
    E = 21000000.0
    v = 0.2
    h = 20
    b = 20
    L = 20

    gamma = 2.54

    coords = [[0.0, 0.0], [L, 0.0], [L, h], [0.0, h]]
    params = Delaunay._strdelaunay(True, True, 1, 30, 1)
    geo = Delaunay(coords, params, 3)
    geo.show()
    plt.show()

    geometria = geo.extrude(b, 5, fast=True)

    O = Elasticity(geometria, E, v, gamma/9.81,
                   verbose=True, solver=LinealEigen)
    O.solve()
    O.exportJSON("Examples/Mesh_tests/Example48.json")
    # a = 0
