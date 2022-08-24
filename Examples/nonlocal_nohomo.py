if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    from FEM.Elasticity2D import PlaneStressNonLocalSparseNonHomogeneous
    from FEM.Geometry.Geometry import Geometry2D
    from FEM.Utils import enmalladoFernando

    E = 2.1*10**6
    v = 0.2
    u = 0.001
    t = 0.5
    l = 0.1
    z1 = 0.5
    LR = 6*l
    P = 1
    a = 5

    nx = 30
    ny = 30

    coords, dicc = enmalladoFernando(a, a, nx, ny)

    geometria = Geometry2D(dicc, coords, ['C2V']*len(dicc), 2, fast=True)

    geometria.generateRegionFromCoords([0, 0], [a, 0])
    geometria.generateRegionFromCoords([a, 0], [a, a])
    geometria.generateRegionFromCoords([a, a], [0, a])
    geometria.generateRegionFromCoords([0, a], [0, 0])
    cb = geometria.cbFromRegion(3, 0, 1)
    cb += geometria.cbFromRegion(3, 0, 2)
    cb += geometria.cbFromRegion(1, u, 1)
    geometria.setCbe(cb)

    # geometria.show(draw_bc=True, label_bc=True)
    # plt.show()

    def af(rho):
        return (0.5/np.pi/l/l/t)*np.exp(-rho)

    O = PlaneStressNonLocalSparseNonHomogeneous(
        geometria, E, v, t, l, z1, Lr=6*l, af=af, verbose=True)
    O.solve()
    O.exportJSON('NonLocalNonHomogeneous.json')
    plt.show()
