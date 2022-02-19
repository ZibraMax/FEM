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
    'Examples/Mesh_tests/Beam_serendipity.json', fast=True)
cbe = geometria.cbFromRegion(1, 0, 1)
cbe += geometria.cbFromRegion(1, 0, 2)
cbe += geometria.cbFromRegion(3, 0, 1)
cbe += geometria.cbFromRegion(3, 0, 2)
geometria.cbe = cbe
O = PlaneStressSparse(geometria, E, v, b, fy=lambda x: -gamma*b, verbose=True)
O.geometry.mask = None
O.solve()
plt.show()
