if __name__ == '__main__':
    from FEM.Geometry import Lineal
    from FEM.Heat1D import Heat1DTransient
    import matplotlib.pyplot as plt
    L = 1
    A = 1
    P = 1
    k = 1
    beta = 0
    Ta = 0
    T0 = 0

    n = 1
    geometria = Lineal(L, n, 1)
    geometria.cbe = [[0, T0]]
    O = Heat1DTransient(geometria, A, P, k, beta, Ta, verbose=False)
    O.set_initial_condition(1.0)
    O.set_alpha(0.0)
    O.solve(t0=0, tf=6, steps=50, dt=0.05, plot=False)
    X, U1L0 = O.postProcess(node=-1)

    L = 1
    A = 1
    P = 1
    k = 1
    beta = 0
    Ta = 0
    T0 = 0

    n = 1
    geometria = Lineal(L, n, 1)
    geometria.cbe = [[0, T0]]
    O = Heat1DTransient(geometria, A, P, k, beta, Ta, verbose=False)
    O.set_initial_condition(1.0)
    O.set_alpha(1.0)
    O.solve(t0=0, tf=6, steps=50, dt=0.05, plot=False)
    X, U1L1 = O.postProcess(node=-1, ax=plt.gca())
    plt.gca().get_lines()[-1].set_color("yellow")

    L = 1
    A = 1
    P = 1
    k = 1
    beta = 0
    Ta = 0
    T0 = 0

    n = 1
    geometria = Lineal(L, n, 1)
    geometria.cbe = [[0, T0]]
    O = Heat1DTransient(geometria, A, P, k, beta, Ta, verbose=False)
    O.set_initial_condition(1.0)
    O.set_alpha(0.5)
    O.solve(t0=0, tf=6, steps=50, dt=0.05, plot=False)
    X, U1L0_5 = O.postProcess(node=-1, ax=plt.gca())
    plt.gca().get_lines()[-1].set_color("red")

    L = 1
    A = 1
    P = 1
    k = 1
    beta = 0
    Ta = 0
    T0 = 0

    n = 2
    geometria = Lineal(L, n, 1)
    geometria.cbe = [[0, T0]]
    O = Heat1DTransient(geometria, A, P, k, beta, Ta, verbose=False)
    O.set_initial_condition(1.0)
    O.set_alpha(0.5)
    O.solve(t0=0, tf=6, steps=50, dt=0.05, plot=False)
    X, U2L0_5 = O.postProcess(node=-1, ax=plt.gca())
    plt.gca().get_lines()[-1].set_color("blue")
    # O.exportJSON("Examples/Mesh_tests/Example51.json")

    L = 1
    A = 1
    P = 1
    k = 1
    beta = 0
    Ta = 0
    T0 = 0

    n = 4
    geometria = Lineal(L, n, 1)
    geometria.cbe = [[0, T0]]
    O = Heat1DTransient(geometria, A, P, k, beta, Ta, verbose=False)
    O.set_initial_condition(1.0)
    O.set_alpha(0.5)
    O.solve(t0=0, tf=6, steps=50, dt=0.05, plot=False)
    X, U4L0_5 = O.postProcess(node=-1, ax=plt.gca())

    L = 1
    A = 1
    P = 1
    k = 1
    beta = 0
    Ta = 0
    T0 = 0

    n = 8
    geometria = Lineal(L, n, 1)
    geometria.cbe = [[0, T0]]
    O = Heat1DTransient(geometria, A, P, k, beta, Ta, verbose=False)
    O.set_initial_condition(1.0)
    O.set_alpha(0.5)
    O.solve(t0=0, tf=6, steps=50, dt=0.05, plot=False)
    X, U8L0_5 = O.postProcess(node=-1, ax=plt.gca())

    L = 1
    A = 1
    P = 1
    k = 1
    beta = 0
    Ta = 0
    T0 = 0

    n = 1
    geometria = Lineal(L, n, 2)
    geometria.cbe = [[0, T0]]
    O = Heat1DTransient(geometria, A, P, k, beta, Ta, verbose=False)
    O.set_initial_condition(1.0)
    O.set_alpha(0.5)
    O.solve(t0=0, tf=6, steps=50, dt=0.05, plot=False)
    X, U1Q0_5 = O.postProcess(node=-1, ax=plt.gca())

    L = 1
    A = 1
    P = 1
    k = 1
    beta = 0
    Ta = 0
    T0 = 0

    n = 2
    geometria = Lineal(L, n, 2)
    geometria.cbe = [[0, T0]]
    O = Heat1DTransient(geometria, A, P, k, beta, Ta, verbose=False)
    O.set_initial_condition(1.0)
    O.set_alpha(0.5)
    O.solve(t0=0, tf=6, steps=50, dt=0.05, plot=False)
    X, U2Q0_5 = O.postProcess(node=-1, ax=plt.gca())

    L = 1
    A = 1
    P = 1
    k = 1
    beta = 0
    Ta = 0
    T0 = 0

    n = 4
    geometria = Lineal(L, n, 2)
    geometria.cbe = [[0, T0]]
    O = Heat1DTransient(geometria, A, P, k, beta, Ta, verbose=False)
    O.set_initial_condition(1.0)
    O.set_alpha(0.5)
    O.solve(t0=0, tf=6, steps=50, dt=0.05, plot=False)
    X, U4Q0_5 = O.postProcess(node=-1, ax=plt.gca())

import numpy as np
SOLUTIONS = np.array([U1L0, U1L1, U1L0_5, U2L0_5, U4L0_5,
                     U8L0_5, U1Q0_5, U2Q0_5, U4Q0_5]).T
np.savetxt("temperaturas_transient.csv", SOLUTIONS[0], delimiter=",", fmt="%s")
