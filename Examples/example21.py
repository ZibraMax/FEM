if __name__ == '__main__':
    import numpy as np
    from FEM.Elasticity2D import PlaneStressSparse
    from FEM.Geometry.Region import Region1D
    from FEM.Geometry import Geometry2D

    b = 120
    h = 160
    t = 0.036
    E = 30*10**(6)
    v = 0.25

    gdls = np.array([[0, 0], [b, 0], [0, h], [b, h]])
    elemento1 = [0, 1, 3, 2]
    dicc = [elemento1]
    tipos = ['C1V']
    regiones = [[0, 1], [1, 3], [3, 2], [2, 0]]
    regions = []
    for reg in regiones:
        regions += [Region1D(gdls[np.ix_(reg)])]
    geometria = Geometry2D(dicc, gdls, tipos, nvn=2,
                           regions=regions, fast=True)
    cbe = geometria.cbFromRegion(3, 0, 1)
    cbe += geometria.cbFromRegion(3, 0, 2)
    geometria.cbe = cbe
    geometria.loadOnRegion(1, fx=lambda s: 10, fy=lambda s: 0)

    O = PlaneStressSparse(geometria, E, v, t)
    O.elementMatrices()
    O.ensembling()
    O.borderConditions()
    O.solveES()
    O.postProcess()
    # plt.show()

    print(O.giveStressPoint(np.array([[60], [80]])))
