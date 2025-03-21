if __name__ == '__main__':
    from FEM.EulerBernoulliBeam import EulerBernoulliBeamNonLinear
    from FEM.Geometry import Lineal
    import matplotlib.pyplot as plt

    L = 100
    geometria = Lineal(L, 100, 1, 3)
    geometria.cbe = [[0, 0.0], [1, 0.0], [
        2, 0.0], [-1, 0.0], [-2, 0.0], [-3, 0.0]]

    W = 10.0
    Fx = 0.0
    I = 1/12
    A = 1
    E = 30000000.0
    EI = E*I
    EA = E*A

    O = EulerBernoulliBeamNonLinear(geometria, EI, EA, Fx, W)
    O.solve()
    plt.close("all")
    for i in range(len(O.solver.solutions)):
        O.solver.setSolution(i, True)
        O.postProcess()
    O.exportJSON("Examples/Mesh_tests/Example41.json")
    plt.show()
    a = 0
