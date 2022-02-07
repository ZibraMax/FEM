import numpy as np
import matplotlib.pyplot as plt
from FEM.Torsion2D import Torsion2D
from FEM.Geometry.Geometry import Geometry

G = 1
phi = 1

geometria = Geometry.loadmsh('Mesh_tests/Web_test.msh')
O = Torsion2D(geometria, G, phi)
O.solve()
plt.show()
