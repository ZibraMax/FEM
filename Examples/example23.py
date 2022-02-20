if __name__ == '__main__':
    from FEM.Geometry import Geometry2D
    from FEM.Elasticity2D import PlaneStrainSparse
    from FEM.Geometry.Region import Region1D
    import matplotlib.pyplot as plt
    import numpy as np
    # P 11.1
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
        print(U)

        def f(x):
            suma = 0
            for i in range(n):
                suma += U[i, 0]*x**i
            return suma
        return f

    parabola = darPolinomio(np.array([0, 10, 20]), np.array(
        [0, b-ancho_en_h10_in/2, b-ancho_en_h20_in/2]))
    coordendas = np.array([
        [0, 0],
        [b/2, 0],
        [b, 0],
        [3*b/2, 0],
        [2*b, 0],

        [parabola(he/2), he/2],
        [b, he/2],
        [2*b-parabola(he/2), he/2],

        [parabola(he), he],
        [b/2+parabola(he)/2, he],
        [b, he],
        [3*b/2-parabola(he)/2, he],
        [2*b-parabola(he), he],

        [parabola(3*he/2), 3*he/2],
        [b, 3*he/2],
        [2*b-parabola(3*he/2), 3*he/2],

        [parabola(2*he), 2*he],
        [b/2+parabola(2*he)/2, 2*he],
        [b, 2*he],
        [3*b/2-parabola(2*he)/2, 2*he],
        [2*b-parabola(2*he), 2*he],

        [parabola(5*he/2), 5*he/2],
        [b, 5*he/2],
        [2*b-parabola(5*he/2), 5*he/2],

        [parabola(3*he), 3*he],
        [b/2+parabola(3*he)/2, 3*he],
        [b, 3*he],
        [3*b/2-parabola(3*he)/2, 3*he],
        [2*b-parabola(3*he), 3*he],

        [parabola(7*he/2), 7*he/2],
        [b, 7*he/2],
        [2*b-parabola(7*he/2), 7*he/2],

        [parabola(4*he), 4*he],
        [b/2+parabola(4*he)/2, 4*he],
        [b, 4*he],
        [3*b/2-parabola(4*he)/2, 4*he],
        [2*b-parabola(4*he), 4*he]])

    elementos = [
        [1, 3, 11, 9, 2, 7, 10, 6],
        [3, 5, 13, 11, 4, 8, 12, 7],
        [9, 11, 19, 17, 10, 15, 18, 14],
        [11, 13, 21, 19, 12, 16, 20, 15],
        [17, 19, 27, 25, 18, 23, 26, 22],
        [19, 21, 29, 27, 20, 24, 28, 23],
        [25, 27, 35, 33, 26, 31, 34, 30],
        [27, 29, 37, 35, 28, 32, 36, 31]]
    tipos = ['C2V']*len(elementos)
    regiones = [[0, 4], [32, 36]]
    regions = []
    for reg in regiones:
        regions += [Region1D(coordendas[np.ix_(reg)])]
    nvn = 2
    geometria = Geometry2D((np.array(elementos)-1).astype(int).tolist(),
                           coordendas, tipos, nvn, regions, fast=True)
    geometria.cbe = geometria.cbFromRegion(0, 0, 1)
    geometria.cbe += geometria.cbFromRegion(0, 0, 2)
    geometria.cbn = [[50, ppx], [51, ppy]]
    geometria.loadOnRegion(1, fy=lambda s: -p0)
    O = PlaneStrainSparse(geometria, E, v, verbose=True)
    O.elementMatrices()
    O.ensembling()
    O.borderConditions()
    O.solveES()
    O.postProcess()

    plt.show()
