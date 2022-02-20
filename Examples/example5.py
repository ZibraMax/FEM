if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from FEM.Torsion2D import Torsion2D
    from FEM.Geometry import Delaunay

    a = 1
    b = 1
    E = 200000
    v = 0.2
    G = E / (2 * (1 + v))
    phi = 1
    vertices = [[0, 0], [a, 0], [a, a], [0, a]]
    params = Delaunay._strdelaunay(
        constrained=True, delaunay=True, a='0.0005', o=2)
    geometria = Delaunay(vertices, params)
    geometria.exportJSON('Examples/Mesh_tests/Square_torsion.json')
    O = Torsion2D(geometria, G, phi)
    O.solve()
    plt.show()
