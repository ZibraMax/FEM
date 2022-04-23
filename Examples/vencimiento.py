if __name__ == '__main__':
    import types
    import numpy as np
    import matplotlib.pyplot as plt
    from FEM.Geometry import Geometry2D
    from FEM.Elasticity2D import PlaneStress
    import triangle as tr

    filename = 'Datos.csv'
    Rawdata = np.genfromtxt(filename, delimiter=',')

    x_pos = Rawdata[:, 0]
    y_pos = Rawdata[:, 2]
    x_dir1 = Rawdata[:, 3]
    y_dir1 = Rawdata[:, 4]*0.5

    coords = np.array([x_pos, y_pos]).T

    original = dict(vertices=coords)
    triangular = tr.triangulate(original, "a10")
    coords = triangular['vertices']
    eles = triangular['triangles']

    geo = Geometry2D(eles.tolist(), coords, types=['T1V']*len(eles), nvn=2)
    # geo.exportJSON('Poaso.json')
    U = []  # np.array([x_dir1, y_dir1])
    for i in range(len(x_dir1)):
        U += [x_dir1[i]]
        U += [y_dir1[i]]
    U = np.array(U)

    E = 97  # MPa
    v = 0.3
    t = 0.65

    O = PlaneStress(geo, E, v, t)
    O.solveFromArray(U, mult=20)
    O.solver.solutions = [U]
    O.exportJSON('Poaso.json')
    # geo.show()
    # plt.plot(*coords.T, 'o')
    plt.show()

    # color1 = Rawdata[:,6]
    # x_dir2 = Rawdata[:,7]
    # y_dir2 = Rawdata[:,8]*0.5
    # color2 = Rawdata[:,10]
