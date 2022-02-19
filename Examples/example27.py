if __name__ == '__main__':
    from FEM.Geometry import Delaunay
    from FEM.Elasticity2D import PlaneStressSparse
    from FEM.Utils.polygonal import giveCoordsCircle
    import matplotlib.pyplot as plt

    coords = [[0, 0], [4, 0], [4, 0.5], [6, 0.5],
              [6, 2.5], [4, 2.5], [4, 3], [0, 3]]
    fillets = [{'start_region': 1, 'end_region': 2, 'r': 0.48, 'n': 5},
               {'start_region': 4, 'end_region': 5, 'r': 0.48, 'n': 5}]
    holes = []
    radi = 0.5
    cent = [2, 1.5]
    vert, seg = giveCoordsCircle(cent, radi, n=50)
    hole = {'center': cent, 'regions': seg, 'vertices': vert}
    holes += [hole]
    params = Delaunay._strdelaunay(
        constrained=True, delaunay=True, a='0.01', o=2)
    geometria = Delaunay(coords, params, nvn=2, fillets=fillets,
                         holes_dict=holes, fast=True)

    geometria.cbe = geometria.cbFromRegion(7, 0, 1)
    geometria.cbe += geometria.cbFromRegion(7, 0, 2)

    geometria.show()
    geometria.exportJSON('Examples/Mesh_tests/pieza_acero.json')
    plt.show()

    E = 29000000
    v = 0.26
    t = 0.5
    p0 = 5000
    p = p0/2
    geometria.loadOnRegion(3, lambda s: p)
    O = PlaneStressSparse(geometria, E, v, t, verbose=True)
    O.solve()
    plt.show()
