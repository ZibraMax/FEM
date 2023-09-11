if __name__ == '__main__':
    from FEM import EDO1D
    from FEM.Geometry import Lineal
    import matplotlib.pyplot as plt
    import numpy as np

    L = 9
    I1 = 0.000276041666666667
    I2 = 0.000341145833333333
    E = 2.10E+08
    W = 100
    nel = 21
    geo = Lineal(L, nel, 2, 1)

    def EI(x):
        EE = (x <= 1/3*L)*E*I1
        EE += (x >= 2/3*L)*E*I1
        EE += (x < 2/3*L)*(x > 1/3*L)*E*I2
        return EE

    def M(x):
        m = W/2*(L*x-x**2)
        return -m
    CBE = [[0, 0],
           [-1, 0]]
    geo.setCbe(CBE)
    def c(x): return 0
    O = EDO1D(geo, EI, c, M, verbose=True)
    O.solve(plot=True)
    np.savetxt("gato.csv", O.K, delimiter=",", fmt="%s")
    np.savetxt("gatoF.csv", O.F, delimiter=",", fmt="%s")
    np.savetxt("gatoUWU.csv", O.U, delimiter=",", fmt="%s")
    O.exportJSON("Ejemplogato.json")
    plt.show()
