import numpy as np
import matplotlib.pyplot as plt
import FEM
from FEM import Mesh

geometria = Mesh.Geometry.loadmsh('Web_test.msh')
geometria.show()
plt.show()
O = FEM.Torsion2D(geometria,1,1)
O.solve()
plt.show()
