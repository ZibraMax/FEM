import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity2D import PlaneStressSparse
from FEM.Geometry import Geometry2D

E = 21000000.0  # MPa
v = 0.2  # m
h = 0.6  # m
b = 0.3  # m
L = 2.5  # m
a = h**2/100
gamma = 23.54
geometria = Geometry2D.importJSON(
    'Examples/Mesh_tests/beam_serendipity.json', fast=True)
geometria.cbe = geometria.cbFromRegion(3, 0, 1)
geometria.cbe += geometria.cbFromRegion(3, 0, 2)
O = PlaneStressSparse(geometria, E, v, b, fy=lambda x: -gamma, verbose=True)
O.solve(plot=False)
O.profile([0.0, 0], [0.0, h])
# O.profile([0,h],[L,h])
plt.show()
