if __name__ == '__main__':
    import json
    import numpy as np
    from FEM.Elasticity3D import NonLocalElasticityFromTensor
    from FEM.Geometry.Geometry import Geometry3D
    from FEM.Solvers import LinealEigen

    c11 = 166.0  # MPa
    c12 = 63.9  # MPa
    c44 = 79.6  # MPa
    C = np.array([
        [c11, c12, c12, 0, 0, 0],
        [c12, c11, c12, 0, 0, 0],
        [c12, c12, c11, 0, 0, 0],
        [0, 0, 0, c44, 0, 0],
        [0, 0, 0, 0, c44, 0],
        [0, 0, 0, 0, 0, c44]])*10**9  # Pa-3

    h = 20.4356  # Armstrong
    b = 20.4356  # Armstrong
    L = 20.4356  # Armstrong

    rho = 2.329  # g/cm³

    l = 0.535
    z1 = 0.5
    Lr = 6*l

    def af(rho):
        return (1/(8*np.pi*l**3))*np.exp(-rho)  # No referencia, sacada a mano

    geometria = Geometry3D.importJSON('piramid.json', fast=True)

    O = NonLocalElasticityFromTensor(
        geometria, C, rho, l, z1, Lr, af, verbose=True, solver=LinealEigen)
    O.solve(path='tetra.csv')
    y = O.geometry.exportJSON()
    pjson = json.loads(y)
    pjson["disp_field"] = O.eigvec.real.T.tolist()
    y = json.dumps(pjson)
    with open("../FEM C++/docs/resources/TETRA2.json", "w") as f:
        f.write(y)
