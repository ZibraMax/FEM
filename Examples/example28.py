if __name__ == '__main__':
    from FEM.Geometry import Geometry2D
    from FEM.Elasticity2D import PlaneStressSparse
    import matplotlib.pyplot as plt

    geometria = Geometry2D.importJSON(
        'Examples/Mesh_tests/Example27.json', fast=True)
    geometria.show()
    plt.show()
    E = 29000000
    v = 0.26
    t = 0.5
    p0 = 5000
    p = p0/2
    geometria.loadOnRegion(3, lambda s: p)
    O = PlaneStressSparse(geometria, E, v, t)
    O.solve()
    O.exportJSON("Examples/Mesh_tests/Example28.json")
    plt.show()
