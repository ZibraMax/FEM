from FEM.Mesh.Geometry import Geometry
from FEM.PlaneStress import PlaneStress
import matplotlib.pyplot as plt
import numpy as np

geometria = Geometry.loadmsh('Mesh_tests/pieza_acero.msh')
geometria.show()
plt.show()
E = 29000000
v = 0.26
t = 0.5
p0 = 5000
p = p0/2
geometria.loadOnSegment(3, lambda s: p)
O = PlaneStress(geometria, E, v, t)
O.solve()
plt.show()
