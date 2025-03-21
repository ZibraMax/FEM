if __name__ == '__main__':
    from FEM import BarAndHingeLinear
    from FEM.Geometry import Geometry3D
    import matplotlib.pyplot as plt
    import numpy as np

    def R(theta):
        return theta*np.pi/180

    def D(theta):
        return theta*180/np.pi

    def kf(theta):
        return 1

    L = 1
    EA = 1

    x_coord = np.sin(np.pi/3)*L
    theta_0 = R(210)
    phi = 2*np.pi-theta_0
    z_coord = np.sin(phi)*x_coord
    x_coord2 = np.cos(phi)*x_coord

    coords = ([[0, -0.5*L, 0],
               [0, 0.5*L, 0],
               [x_coord, 0, 0],
               [x_coord2, 0, z_coord]])

    elements = [[0, 1],
                [1, 2],
                [2, 0],
                [0, 3],
                [1, 3],
                [3, 0, 1, 2]]

    types = ['L1V']*(len(elements)-1) + ['OH']*1

    geo = Geometry3D(elements, coords, types, 3, fast=True)
    for node in [0, 1, 2]:
        geo.cbe += [[node*3, 0], [node*3+1, 0], [node*3+2, 0]]
    O = BarAndHingeLinear(geo, EA, 1, kf)
    O.addLoadNode(3, [0.0, 0, 1.0])
    O.solve()
    O.exportJSON('./Examples/Mesh_tests/Bar_and_hinge.json')
    a = 0
