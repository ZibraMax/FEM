if __name__ == '__main__':
    import scipy
    import json
    import numpy as np
    from FEM.Elasticity3D import Elasticity
    from FEM.Geometry.Geometry import Geometry3D

    E = 21000000.0
    v = 0.2
    h = 1
    b = 1
    L = 1

    gamma = 23.54

    geometria = Geometry3D.importJSON('piramid.json', fast=True)

    O = Elasticity(geometria, E, v, gamma, verbose=True)
    O.solve()
    print('sssssssss')
    O.M = O.M.tocsc()
    print('sssssssss')
    O.K = O.K.tocsc()
    print('sssssssss')
    omh = scipy.sparse.linalg.spsolve(O.M, O.K)
    print('sssssssss')
    eiv, eigvec = scipy.sparse.linalg.eigs(omh, k=30, which='SR')
    print('Open c++')
    y = O.geometry.exportJSON()
    pjson = json.loads(y)
    pjson["disp_field"] = eigvec.real.T.tolist()
    y = json.dumps(pjson)
    with open("../FEM C++/docs/resources/exported2.json", "w") as f:
        f.write(y)
    np.savetxt('a.csv', O.U, delimiter=',', fmt='%s')
    a = 0
