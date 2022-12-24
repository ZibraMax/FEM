if __name__ == '__main__':
    from FEM.Geometry import Delaunay
    from FEM.Elasticity2D import PlaneStrainSparse
    import matplotlib.pyplot as plt
    import numpy as np

    E = 30*10**(5)
    v = 0.25

    b = 10
    h = 20
    he = h/4
    ancho_en_h10_in = 18
    ancho_en_h20_in = 10
    p0 = 200
    pp = 1000
    ppx = pp*3/5
    ppy = -pp*4/5

    def darPolinomio(X, Y):
        n = len(X)
        A = np.zeros([n, n])
        B = np.zeros([n, 1])
        for i in range(n):
            for j in range(n):
                A[i, j] = X[i]**j
            B[i, 0] = Y[i]

        U = np.linalg.solve(A, B)

        def f(x):
            suma = 0
            for i in range(n):
                suma += U[i, 0]*x**i
            return suma
        return f

    n = 10

    parabola = darPolinomio(np.array([0, 10, 20]), np.array(
        [0, b-ancho_en_h10_in/2, b-ancho_en_h20_in/2]))

    c = [
        [0, 0],
        [2*b, 0]]

    for i in range(n):
        x = 2*b-parabola(h/n*i)
        y = h/n*i
        c += [[x, y]]

    c += [[2*b-parabola(4*he), 4*he],
          [parabola(4*he), 4*he]]
    for i in reversed(range(n)):
        x = parabola(h/n*i)
        y = h/n*i
        c += [[x, y]]
    params = Delaunay._strdelaunay(
        constrained=True, delaunay=True, a='0.5', o=2)
    geometria = Delaunay(c, params, nvn=2, fast=True)
    geometria.generateRegionFromCoords([0, 0], [2*b, 0])
    geometria.generateRegionFromCoords(
        [2*b-parabola(4*he), 4*he], [parabola(4*he), 4*he])
    geometria.cbe = geometria.cbFromRegion(-2, 0, 1)
    geometria.cbe += geometria.cbFromRegion(-2, 0, 2)
    cbn = geometria.generateBCFromCoords(parabola(h/2), h/2, ppx, nv=1)
    cbn += geometria.generateBCFromCoords(parabola(h/2), h/2, ppx, nv=2)
    geometria.cbn = cbn
    geometria.loadOnRegion(-1, fy=lambda s: -p0)
    O = PlaneStrainSparse(geometria, E, v, verbose=True)
    O.elementMatrices()
    O.ensembling()
    O.borderConditions()
    O.solveES()
    O.postProcess()
    O.exportJSON("Examples/Mesh_tests/Example24.json")

    plt.show()
