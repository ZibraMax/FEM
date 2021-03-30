import numpy as np
import matplotlib.pyplot as plt
from FEM.PlaneStress import PlaneStress
from FEM.Mesh.Geometry import Geometry

E = 21000000.0  # MPa
v = 0.2  # m
h = 0.6  # m
b = 0.3  # m
L = 2.5  # m
a = h**2/100
gamma = 23.54
geometria = Geometry.loadmsh('Mesh_tests/Beam_serendipity.msh')
cbe = geometria.cbFromSegment(1, 0, 1)
cbe += geometria.cbFromSegment(1, 0, 2)
cbe += geometria.cbFromSegment(3, 0, 1)
cbe += geometria.cbFromSegment(3, 0, 2)
geometria.cbe = cbe
geometria.loadOnSegment(2, fy=lambda s: -23.54*b*h)
O = PlaneStress(geometria, E, v, b)
O.geometry.show()
plt.show()
O.solve()
n = len(O.U)
pares = np.linspace(0, n-1, n)
print(np.max(np.abs(O.U[pares % 2 == 0])), np.max(np.abs(O.U[pares % 2 == 1])))
plt.show()
