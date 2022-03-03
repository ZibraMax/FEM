if __name__ == '__main__':
    import json
    import numpy as np
    from FEM.Elasticity3D import ElasticityFromTensor
    from FEM.Geometry.Geometry import Geometry3D
    from FEM.Solvers import LinealEigen

    c11 = 223.1  # MPa
    c12 = 63.9  # MPa
    c44 = 79.6  # MPa
    C = np.array([
        [c11, c12, c12, 0, 0, 0],
        [c12, c11, c12, 0, 0, 0],
        [c12, c12, c11, 0, 0, 0],
        [0, 0, 0, c44, 0, 0],
        [0, 0, 0, 0, c44, 0],
        [0, 0, 0, 0, 0, c44]])*10**9  # Pa-3

    rho = 2.329  # g/cm³

    # geometria = Geometry3D.importJSON(
    #     '../FEM C++/docs/resources/SPHERE_FERNANDO.json', fast=True)
    geometria = Geometry3D.importJSON(
        'SPHERE_FERNANDO.json', fast=True)

    O = ElasticityFromTensor(
        geometria, C, rho, verbose=True, solver=LinealEigen)
    O.solve(path='sphere_si_local.csv')
    O.exportJSON('../FEM C++/docs/resources/SPHERE202020.json')
