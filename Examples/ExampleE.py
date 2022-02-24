

from operator import ge
from FEM.Geometry.Region import Region2D


if __name__ == '__main__':
    import scipy
    import json
    import numpy as np
    import matplotlib.pyplot as plt
    from FEM.Elasticity3D import Elasticity
    from FEM.Geometry import Geometry3D, Delaunay

    E = 21000000.0
    v = 0.2
    h = 0.6
    b = 0.3
    L = 2.5

    gamma = 23.54

    coords = [[0.0, 0.0], [L, 0.0], [L, h], [0.0, h]]
    params = Delaunay._strdelaunay(True, True, 0.01, 30, 1)
    geo = Delaunay(coords, params, 1)
    geo.show()
    plt.show()

    base = Region2D(
        np.array([[0.0, 0.0, 0.0], [0.0, h, 0.0], [0.0, h, b], [0.0, 0.0, b]]))
    geometria = geo.extrude(b, 5, fast=True, regions=[base])
    cb = geometria.cbFromRegion(0, 0.0, 1)
    cb += geometria.cbFromRegion(0, 0.0, 2)
    cb += geometria.cbFromRegion(0, 0.0, 3)
    geometria.setCbe(cb)

    O = Elasticity(geometria, E, v, gamma/9.81,
                   fy=lambda x: -gamma, verbose=True)
    O.solve()
    print('A')
    omh = scipy.sparse.linalg.spsolve(O.M, O.K)
    print('A')
    eiv, eigvec = scipy.sparse.linalg.eigs(omh, k=30, which='SR')
    print('Open c++')
    y = O.geometry.exportJSON()
    pjson = json.loads(y)
    pjson["disp_field"] = eigvec.real.T.tolist()
    y = json.dumps(pjson)
    with open("../FEM C++/docs/resources/exported3.json", "w") as f:
        f.write(y)
    # np.savetxt('a.csv', O.U, delimiter=',', fmt='%s')
    # a = 0
