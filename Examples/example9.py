import numpy as np
import matplotlib.pyplot as plt
from FEM.Mesh.Geometry import Geometry
from FEM.Elasticity2D import PlaneStress

E = 21000000.0  # MPa
v = 0.2  # m
h = 0.6  # m
b = 0.3  # m
L = 2.5  # m
a = h**2/100
gamma = 23.54
geometria = Geometry.loadmsh('Mesh_tests/Beam_serendipity.msh')
cbe = geometria.generateBCFromCoords(0, h/2, 0, 1)
cbe += geometria.generateBCFromCoords(0, h/2, 0, 2)
cbe += geometria.generateBCFromCoords(L, h/2, 0, 1)
cbe += geometria.generateBCFromCoords(L, h/2, 0, 2)
geometria.cbe = cbe
O = PlaneStress(geometria, E, v, b, fy=lambda x: -gamma)
O.solve()
plt.show()
