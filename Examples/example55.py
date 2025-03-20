if __name__ == '__main__':
    from FEM import TrussLinear, LinealElement
    from FEM.Geometry import Geometry3D
    import matplotlib.pyplot as plt
    import numpy as np

    LX = 1
    LY = 1
    LZ = 0.5
    EA = 1

    coords = [[0, -0.5*LY, 0],
              [0, 0.5*LY, 0],
              [LX, 0, 0],
              [-0.5*LX, 0, LZ]]

    conectivity = [[0, 1],
                   [1, 2],
                   [2, 0],
                   [0, 3],
                   [1, 3],
                   [3, 2]]

    types = ['L1V']*len(conectivity)

    geo = Geometry3D(conectivity, coords, types, 3, fast=True)
    for node in [0, 1, 2]:
        geo.cbe += [[node*3, 0], [node*3+1, 0], [node*3+2, 0]]
    O = TrussLinear(geo, EA, 1)
    O.addLoadNode(3, [1.0, 0, 0])
    O.solve()
    O.exportJSON('truss.json')
    a = 0
