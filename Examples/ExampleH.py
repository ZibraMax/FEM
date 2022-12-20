if __name__ == '__main__':
    import json
    import numpy as np
    import matplotlib.pyplot as plt
    from FEM.Elasticity3D import Elasticity
    from FEM.Geometry.Geometry import Geometry3D
    FILENAME = "Examples\Mesh_tests\SPHERE_FERNANDO.json"
    E = 21000000.0
    v = 0.2
    h = 0.6
    b = 0.3
    L = 3.5
    gamma = 23.54

    geometria = Geometry3D.importJSON(FILENAME, fast=True)
    e = geometria.elements[549]
    geometria.OctTree.query_range_point_radius(
        e.T(e.center.T)[0].flatten(), 11, plot=True)
