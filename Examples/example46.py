if __name__ == '__main__':
    import json
    import numpy as np
    import matplotlib.pyplot as plt
    from FEM.Elasticity3D import Elasticity
    from FEM.Geometry.Geometry import Geometry3D

    E = 21000000.0
    v = 0.2
    h = 0.6
    b = 0.3
    L = 3.5
    gamma = 23.54

    _a = L
    _b = h
    _c = b

    nx = 100
    ny = 10
    nz = 10

    dx = _a/nx
    dy = _b/ny
    dz = _c/nz

    coords = []

    for i in range(nx+1):
        x = i*dx
        for j in range(ny+1):
            y = j*dy
            for k in range(nz+1):
                z = k*dz
                coords += [[x, y, z]]

    dicc = []

    def node(i, j, k): return i*(ny+1)*(nz+1)+j*(nz+1)+k

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                node1 = node(i, j, k)
                node2 = node(i+1, j, k)
                node3 = node(i+1, j+1, k)
                node4 = node(i, j+1, k)

                node5 = node(i, j, k+1)
                node6 = node(i+1, j, k+1)
                node7 = node(i+1, j+1, k+1)
                node8 = node(i, j+1, k+1)

                dicc += [[node1, node2, node3, node4,
                          node5, node6, node7, node8]]

    def fy(x): return -gamma

    geometria = Geometry3D(dicc, coords, ["B1V"]*len(dicc), nvn=3, fast=True)
    cbe = []
    for i in range(len(coords)):
        if 0.0 == coords[i][0]:
            cbe += [[i*3, 0.0]]
            cbe += [[i*3+1, 0.0]]
            cbe += [[i*3+2, 0.0]]
    geometria.cbe = cbe

    O = Elasticity(geometria, E, v, gamma, fy=fy, verbose=True)
    e = O.elements[-50]
    coords = e.coords
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot(*coords.T, '-o', color='black')
    # plt.show()
    O.solve()
    O.exportJSON("Examples/Mesh_tests/Example46.json")
    a = 0
