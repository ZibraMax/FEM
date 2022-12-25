if __name__ == '__main__':
    import FEM
    import matplotlib.pyplot as plt

    b = 30
    h = 30
    k = 311
    kx = 2*k
    ky = k

    params = FEM.Delaunay._strdelaunay(True, True, 2, 30, 2)

    vertices = [[0.0, 0.0], [b, 0.0], [b, 20],
                [20, b], [0, b]]

    geo = FEM.Delaunay(vertices=vertices, params=params, nvn=1)

    cb = geo.cbFromRegion(0, 0)
    cb += geo.cbFromRegion(1, 100)
    cb += geo.cbFromRegion(2, 100)
    cb += geo.cbFromRegion(3, 100)

    geo.setCbe(cb)

    o = FEM.Heat2D(geo, kx, ky, verbose=True)
    o.solve()
    geo.show()
    plt.show()
    o.exportJSON("Examples/Mesh_tests/Example50.json")
