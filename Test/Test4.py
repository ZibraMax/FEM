import numpy as np
import matplotlib.pyplot as plt
import FEM
from FEM import Mesh

G = 1
phi =1

geometria = Mesh.Geometry.loadmsh('Mesh_tests/Web_test.msh')
O = FEM.Torsion2D(geometria,G,phi)
O.solve()
plt.show()
