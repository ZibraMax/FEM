if __name__ == '__main__':
    from FEM.Heat2D import Heat2D
    from FEM.Geometry import Delaunay
    from FEM.Utils.polygonal import giveCoordsCircle
    import matplotlib.pyplot as plt

    coords = [[0, 0], [4, 0], [4, 0.5], [6, 0.5],
              [6, 2.5], [4, 2.5], [4, 3], [0, 3]]
    fillets = [{'start_region': 4, 'end_region': 5, 'r': 0.48, 'n': 5}, {'start_region': 1, 'end_region': 2, 'r': 0.48, 'n': 5},
               ]
    holes = []
    radi = 0.5
    cent = [2, 1.5]
    vert, seg = giveCoordsCircle(cent, radi, n=70)
    hole = {'center': cent, 'regions': seg, 'vertices': vert}
    holes += [hole]

    radi = 0.5/2
    cent = [5, 1.5]
    vert, seg = giveCoordsCircle(cent, radi, n=40)
    hole = {'center': cent, 'regions': seg, 'vertices': vert}
    holes += [hole]

    params = Delaunay._strdelaunay(
        constrained=True, delaunay=True, a='0.01', o=2)
    geometria = Delaunay(coords, params, nvn=1,
                         fillets=fillets, holes_dict=holes)

    cbe = geometria.cbOnHole(0, 100)
    cbe += geometria.cbOnHole(1, 0)
    geometria.setCbe(cbe)

    kx = 15/2.54
    ky = 15/2.54
    Ta = 20
    beta = 5
    O = Heat2D(geometria, kx, ky)
    O.geometry.show(draw_bc=True, draw_segs=False)
    plt.show()
    O.solve()
    O.exportJSON("pieza_acero_heat2d.json")
    plt.show()
