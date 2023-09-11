if __name__ == '__main__':
    from FEM import EulerBernoulliBeam
    from FEM.Geometry import Lineal
    import matplotlib.pyplot as plt
    import numpy as np

    L = 9
    I1 = 0.000380313
    I2 = 0.000445417
    E = 2.10E+08
    W = -100
    nel = 100
    geo = Lineal(L, nel, 1, 2)

    def EI(x):
        EE = (x <= 1/3*L)*E*I1
        EE += (x >= 2/3*L)*E*I1
        EE += (x < 2/3*L)*(x > 1/3*L)*E*I2
        return EE

    O = EulerBernoulliBeam(geo, EI, 0, W, verbose=True)
    CBE = [[0, 0],
           [-2, 0]]
    O.cbe = CBE
    O.solve(plot=False)
    plt.show()
    U1 = O.U.flatten()[::2]

    geo = Lineal(L, nel, 1, 2)

    O = EulerBernoulliBeam(geo, E*I1, 0, W, verbose=True)
    CBE = [[0, 0],
           [-2, 0]]
    O.cbe = CBE
    O.solve(plot=False)
    U2 = O.U.flatten()[::2]
    X = np.linspace(0, L, len(U1))
    plt.plot(X, U1, label="Distintas inercias")
    plt.plot(X, U2, label="Mismas inercias")
    plt.grid()
    plt.legend()
    plt.show()
