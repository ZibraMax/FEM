if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from FEM.Geometry import Geometry2D
    from FEM.Elasticity2D import PlaneStressSparse
    from FEM.Utils import enmalladoFernando

    E = 21000000.0  # MPa
    v = 0.2  # m
    h = 0.6  # m
    b = 0.3  # m
    L = 2.5  # m
    a = h**2/100
    gamma = 23.54

    coords, dicc = enmalladoFernando(L, h, 130, 30)
    geo = Geometry2D(dicc, coords, ['C2V']*len(dicc), nvn=2, fast=True)
    geo.generateRegionFromCoords([0.0, 0.0], [L, 0.0])
    geo.generateRegionFromCoords([L, 0.0], [L, h])
    geo.generateRegionFromCoords([L, h], [0.0, h])
    geo.generateRegionFromCoords([0.0, h], [0.0, 0.0])
    cbe = geo.generateBCFromCoords(0, h/2, 0, 1)
    cbe += geo.generateBCFromCoords(0, h/2, 0, 2)
    cbe += geo.generateBCFromCoords(L, h/2, 0, 1)
    cbe += geo.generateBCFromCoords(L, h/2, 0, 2)
    geo.setCbe(cbe)
    O = PlaneStressSparse(geo, E, v, b,
                          fy=lambda x: -gamma, verbose=True)
    O.solve()
    O.exportJSON("Examples/Mesh_tests/Example9.json")
    plt.show()
