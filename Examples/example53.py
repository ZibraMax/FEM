if __name__ == '__main__':
    from FEM.Geometry import Lineal
    from FEM.Heat1D import Heat1DTransient
    import matplotlib.pyplot as plt
    L = 0.1
    A = 0.001*0.005
    P = 2*(0.005+0.001)
    k = 385
    beta = 25
    Ta = 20
    T0 = 100
    # TODO hay algo mal

    n = 30
    geometria = Lineal(L, n, 2)
    O = Heat1DTransient(geometria, A, P, k, beta, Ta, verbose=False)
    O.set_initial_condition(0.0)
    O.set_alpha(0.5)
    O.cbe = [[0, T0]]
    O.defineConvectiveBoderConditions(-1)
    O.solve(t0=0, tf=10, steps=400, plot=False)
    O.postProcess()
    plt.show()
    O.exportJSON("Examples/Mesh_tests/Example53.json")
