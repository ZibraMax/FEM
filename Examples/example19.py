if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    from FEM.Elasticity2D import PlaneStressSparse
    from FEM.Geometry import Geometry2D

    E = 21000000.0  # MPa
    v = 0.2  # m
    h = 0.6  # m
    b = 0.3  # m
    L = 2.5  # m
    a = h**2/100
    gamma = 23.54
    geometria = Geometry2D.importJSON(
        'Examples/Mesh_tests/Example9.json', fast=True)
    cbe = geometria.cbFromRegion(1, 0, 1)
    cbe += geometria.cbFromRegion(1, 0, 2)
    cbe += geometria.cbFromRegion(3, 0, 1)
    cbe += geometria.cbFromRegion(3, 0, 2)
    geometria.cbe = cbe
    geometria.loadOnRegion(2, fy=lambda s: -23.54*b*h)
    O = PlaneStressSparse(geometria, E, v, b, verbose=True)
    O.geometry.show()
    O.geometry.mask = None
    plt.show()
    O.solve()
    n = len(O.U)
    pares = np.linspace(0, n-1, n)
    print(np.max(np.abs(O.U[pares % 2 == 0])),
          np.max(np.abs(O.U[pares % 2 == 1])))
    O.exportJSON("Examples/Mesh_tests/Example19.json")
    plt.show()
