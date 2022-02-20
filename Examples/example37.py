if __name__ == '__main__':
    from FEM.Heat2D import Heat2D
    from FEM.Geometry import Delaunay
    from FEM.Utils.polygonal import giveCoordsCircle
    import numpy as np
    import matplotlib.pyplot as plt
    # https://github.com/Samson-Mano/2D_Heat_transfer
    k = 2
    t0 = 140

    beta = 1.5
    Ta = 20

    radi = 4
    cent = [0, 0]
    vert_orig, seg_orig = giveCoordsCircle(cent, radi, n=150)

    vert, seg = giveCoordsCircle([0, 1], 1, n=50)
    holes = [{'center': [0, 1], 'regions': seg, 'vertices': vert}]
    params = Delaunay._strdelaunay(
        constrained=True, delaunay=True, a='0.03', o=1)
    geometria = Delaunay(vert_orig, params, 1, holes_dict=holes)
    geometria.setCbe(geometria.cbOnHole(0, t0, 1))

    nodes = geometria.gdls
    n = len(nodes)
    triangles = geometria.dictionary
    piramides = []
    m = 5
    h = 1.0

    dz = h/(m-1)

    dddnodes = np.zeros([m*n, 3])
    for i in range(m):
        dddnodes[n*(i):n*(i+1), :2] = nodes
        dddnodes[n*(i):n*(i+1), -1] = i*dz

    for i in range(m-1):
        for t in triangles:
            t = np.array(t)
            nodossup = n*(i) + t
            nodosinf = n*(i+1) + t
            p = nodossup.tolist()+nodosinf.tolist()

            piramides += [[p[2], p[5], p[0], p[1]]]
            piramides += [[p[0], p[5], p[3], p[4]]]
            piramides += [[p[0], p[5], p[4], p[1]]]

    O = Heat2D(geometria, k, k)
    for i in range(len(seg_orig)):
        O.defineConvectiveBoderConditions(i, beta, Ta)
    O.geometry.show()
    plt.show()
    O.solve()
    plt.show()
