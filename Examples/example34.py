if __name__ == '__main__':
    from FEM.Geometry import Delaunay
    from FEM.Elasticity2D import PlaneStressSparse
    from FEM.Utils.polygonal import giveCoordsCircle
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import gridspec

    b = 5
    h = 11.5

    dh = 3/4
    t = 7/16

    n = 3
    E = 29000000
    v = 0.32
    vu = 1250
    coords = [[0, 0], [b, 0], [b, h], [0, h]]
    holes = []
    for i in range(n):
        radi = dh/2
        cent = [b/2, h/(n+1)*(i+1)]
        vert, seg = giveCoordsCircle(cent, radi, n=20)
        hole = {'center': cent, 'regions': seg, 'vertices': vert}
        holes += [hole]

    params = Delaunay._strdelaunay(
        constrained=True, delaunay=True, a='0.05', o=2)
    geometria = Delaunay(coords, params, nvn=2, holes_dict=holes, fast=True)

    cbe = geometria.cbFromRegion(3, 0, 1)
    cbe += geometria.cbFromRegion(3, 0, 2)
    geometria.setCbe(cbe)

    for i in range(n):
        geometria.loadOnHole(
            i, 2*np.pi, np.pi, fy=lambda s: -vu/(2*np.pi*dh/2/2))
    O = PlaneStressSparse(geometria, E, v, t, verbose=True)
    O.geometry.show(draw_bc=True)
    plt.show()
    O.solve(plot=False)
    gss = gridspec.GridSpec(1, 4)
    gs = [gss[0, 0], gss[0, 1], gss[0, 2], gss[0, 3]]
    O.postProcess(gs=gs)
    plt.show()
