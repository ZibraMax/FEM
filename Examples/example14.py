if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from FEM.Elasticity2D import PlaneStressSparse
    from FEM.Geometry.Geometry import Geometry2D

    E = 21000000.0  # MPa
    v = 0.2  # m
    h = 0.6  # m
    b = 0.3  # m
    L = 2.5  # m
    a = h**2/100
    gamma = 23.54
    geometria = Geometry2D.importJSON(
        'Examples/Mesh_tests/Example10.json', fast=True)
    geometria.cbe = geometria.cbFromRegion(3, 0, 1)
    geometria.cbe += geometria.cbFromRegion(3, 0, 2)
    O = PlaneStressSparse(geometria, E, v, b,
                          fy=lambda x: -gamma, verbose=True)
    O.solve()
    O.exportJSON("Examples/Mesh_tests/Example14.json")
    plt.show()
