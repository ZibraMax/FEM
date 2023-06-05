if __name__ == '__main__':
    from FEM.Geometry import Lineal
    from FEM.Heat1D import Heat1DTransient

    L = 0.1
    A = 0.001*0.005
    P = 2*(0.005+0.001)
    k = 385
    beta = 25
    Ta = 20
    T0 = 100

    n = 40
    geometria = Lineal(L, n, 1)
    O = Heat1DTransient(geometria, A, P, k, beta, Ta, verbose=True)
    O.cbe = [[0, T0]]
    O.cbe = [[-1, 0]]
    O.solve(t0=0, tf=10, steps=100)
    O.exportJSON("Examples/Mesh_tests/Example51.json")
