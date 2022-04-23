if __name__ == '__main__':
    import numpy as np
    import json
    import tetgen
    import pyvista as pv
    import scipy
    import logging
    from FEM.Geometry import Geometry3D
    from FEM import Elasticity
    pv.set_plot_theme('document')
    sphere = pv.Sphere()
    tet = tetgen.TetGen(sphere)
    coords, dicc = tet.tetrahedralize(order=1, mindihedral=20, minratio=1.5)
    grid = tet.grid
    grid.plot(show_edges=True)
    geo = Geometry3D(dicc.tolist(), coords, [
                     'TE1V']*len(dicc), nvn=3, fast=True)
    geo.exportJSON('SPHERE.json')
    print(len(geo.elements), len(geo.gdls)*3)

    E = 21000000.0
    v = 0.2
    gamma = 23.54

    O = Elasticity(geo, E, v, gamma, verbose=True)
    O.solve()
    omh = scipy.sparse.linalg.spsolve(O.M.tocsc(), O.K.tocsc())
    logging.info('Solved M-1K')
    eiv, eigvec = scipy.sparse.linalg.eigs(O.K, k=10, M=O.M, which='SR')
    logging.info('Eigenvalues found.')

    y = O.geometry.exportJSON()
    pjson = json.loads(y)
    pjson["disp_field"] = eigvec.real.T.tolist()
    y = json.dumps(pjson)
    with open("../FEM C++/docs/resources/SPHERE.json", "w") as f:
        f.write(y)
    logging.info('File exported')

    # ledge = 10.0
    # height = ledge*6**0.5/3

    # v = np.array([[0, 0, 0], [ledge, 0, 0.0],
    #               [0, ledge, 0], [0, 0, height]])
    # f = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    # tgen = tetgen.TetGen(v, f)
    # coords, dicc = tgen.tetrahedralize(switches='-pqa10i')
    # tgen.plot(show_edges=True)
