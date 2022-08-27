if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    from FEM.Elasticity2D import PlaneStressNonLocalSparseNonHomogeneous
    from FEM.Geometry import Geometry2D, Region2D
    from FEM.Utils import enmalladoFernando

    E0 = 2.1*10**6
    E1 = 0.4*E0
    E2 = E0
    v = 0.2
    u = 0.001
    t = 0.5
    l = 0.1
    alpha = 50
    # Para que las integrales de Gamma sean mas exactas es necesario un mayor Lr.
    # Ponerle mucho puede ser problemático porque hace que los gamma sean mayores a 1
    Lr = 6*l
    a = 5
    h = 5

    nx = 30
    ny = 30

    coords, dicc = enmalladoFernando(a, h, nx, ny)

    geometria = Geometry2D(dicc, coords, ['C2V']*len(dicc), 2, fast=True)

    geometria.generateRegionFromCoords([0, 0], [a, 0])
    geometria.generateRegionFromCoords([a, 0], [a, a])
    geometria.generateRegionFromCoords([a, a], [0, a])
    geometria.generateRegionFromCoords([0, a], [0, 0])
    cb = geometria.cbFromRegion(3, 0, 1)
    cb += geometria.cbFromRegion(3, 0, 2)
    cb += geometria.cbFromRegion(1, u, 1)
    geometria.setCbe(cb)

    re1 = Region2D(np.array([[2.0, 2.0], [3.0, 2.0], [3.0, 3.0], [2.0, 3.0]]))
    re2 = Region2D(np.array([[0.0, 0.0], [a, 0.0], [a, h], [0.0, h]]))

    geometria.addRegions([re1, re2])

    E = []
    for i in range(len(geometria.elements)):
        if re1.isBetween(geometria.centroids[i][0]):
            E += [E1]
        else:
            E += [E2]

    # geometria.show(draw_bc=True, label_bc=True)
    # plt.show()

    def af(rho):
        return (0.5/np.pi/l/l/t)*np.exp(-rho)

    O = PlaneStressNonLocalSparseNonHomogeneous(
        geometria, E, v, t, l, alpha, Lr, af, verbose=True)
    O.solve(plot=False)
    O.postProcess(mult=10)
    O.exportJSON('NonLocalNonHomogeneousPlate.json')
    plt.show()
