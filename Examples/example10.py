if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from FEM.Elasticity2D import PlaneStressSparse
    from FEM.Geometry import Geometry2D, Delaunay
    import numpy as np

    E = 21000000.0  # MPa
    v = 0.2  # m
    h = 0.6  # m
    b = 0.3  # m
    L = 2.5  # m
    a = h**2/100
    gamma = 23.54

    coords = np.array([[0.0, 0.0], [L, 0.0], [L, h], [0.0, h]])
    params = Delaunay._strdelaunay(a=0.01, q=30, o=2)
    geo = Delaunay(coords, params, nvn=2, fast=True)
    cb = geo.cbFromRegion(3, 0.0, 1)
    cb += geo.cbFromRegion(3, 0.0, 2)
    geo.setCbe(cb)
    geo.exportJSON('Examples/Mesh_tests/triang_beam.json')
    geo.show()
    plt.show()

    geometria = Geometry2D.importJSON(
        'Examples/Mesh_tests/triang_beam.json', fast=True)
    O = PlaneStressSparse(geometria, E, v, b,
                          fy=lambda x: -gamma, verbose=True)
    O.solve()
    O.exportJSON('TRIANG_BEAM.json')
    plt.show()
