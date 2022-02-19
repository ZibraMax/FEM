if __name__ == '__main__':
    from FEM.Elasticity2D import PlaneStressOrthotropic as pso
    from FEM.Geometry import Geometry2D
    import matplotlib.pyplot as plt
    import numpy as np
    from FEM.Geometry.Region import Region1D

    a = 120.0
    b = 160.0
    p0 = 10.0

    E1 = 31*10**6
    E2 = 2.7*10**6
    G12 = 0.75*10**6
    v12 = 0.28
    t = 0.036

    gdls = np.array([[0.0, 0.0], [a, 0.0], [a, b], [0, b]])
    diccs = [[0, 1, 2], [0, 2, 3]]
    types = ['T1V', 'T1V']
    nvn = 2
    regions = [[0, 1], [1, 2], [2, 3], [3, 0]]
    regs = []
    for r in regions:
        regs += [Region1D(gdls[np.ix_(r)])]

    geometry = Geometry2D(dictionary=diccs, gdls=gdls,
                          types=types, nvn=nvn, regions=regs)

    cb = geometry.cbFromRegion(3, 0.0, 1)
    cb += geometry.cbFromRegion(3, 0.0, 2)
    geometry.setCbe(cb)
    def f(s): return p0

    geometry.loadOnRegion(1, f)
    # geometry.show(draw_labels=True, draw_bc=True, label_bc=True)
    # plt.show()
    OFEM = pso(geometry, E1, E2, G12, v12, t)
    OFEM.solve()
    plt.show()
