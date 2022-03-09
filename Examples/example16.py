if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    from FEM.Elasticity2D import PlaneStressNonLocalSparse
    from FEM.Geometry import Geometry2D

    E = 2.1*10**6
    v = 0.2
    u = 0.001
    a = 5
    t = 0.5
    l = 0.1
    z1 = 0.5
    geometria = Geometry2D.importJSON(
        'Examples/Mesh_tests/rect2.json', fast=True)
    cb = geometria.cbFromRegion(3, 0, 1)
    cb += geometria.cbFromRegion(3, 0, 2)
    cb += geometria.cbFromRegion(1, u, 1)
    geometria.setCbe(cb)
    geometria.show(draw_bc=True, label_bc=True)
    plt.show()

    def af(l0, rho):
        return l0*np.exp(-rho)

    O = PlaneStressNonLocalSparse(
        geometria, E, v, t, l, z1, Lr=6*l, af=af, verbose=True)
    O.solve(path="fast.csv")
    plt.show()

    _X, U1, U2, U3, U = O.profile([0, 0.019], [a, 0.019])
    # np.savetxt('a2gt.csv', [_X, U1], delimiter=',')
    plt.show()
