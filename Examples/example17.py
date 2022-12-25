if __name__ == '__main__':
    from FEM.Geometry import Lineal
    from FEM.Heat1D import Heat1D

    L = 0.1
    A = 0.001*0.005
    P = 2*(0.005+0.001)
    k = 385
    beta = 25
    Ta = 20
    T0 = 100

    n = 4
    geometria = Lineal(L, n, 1)
    O = Heat1D(geometria, A, P, k, beta, Ta)
    O.cbe = [[0, T0]]
    O.defineConvectiveBoderConditions(-1)
    O.solve()
    O.exportJSON("Examples/Mesh_tests/Example17.json")
