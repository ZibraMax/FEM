if __name__ == '__main__':
    import numpy as np
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
