if __name__ == '__main__':
    from FEM import EDO1D
    from FEM.Geometry import Geometry1D
    import matplotlib.pyplot as plt
    import numpy as np
    E = 20000000
    r1 = 0.5
    r2 = 0.3
    A1 = r1**2*np.pi
    A2 = r2**2*np.pi

    def EA(x):
        EE = (x <= 3)*E*A1
        EE += (x > 3)*E*A2
        return EE

    geo = Geometry1D([[0, 1],
                      [1, 2],
                      [2, 3],
                      [3, 4]],

                     [[0.0],
                      [1.5],
                      [3],
                      [4],
                      [5]],

                     ["L1V"]*4, 1)
    geo.setCbe([[0, 0]])
    geo.cbn = [[-1, -850]]
    O = EDO1D(geo, EA, lambda x: 0, lambda x: 0)
    O.solve(plot=True)
    O.exportJSON("Ejemplo_axial.json")
    plt.show()
