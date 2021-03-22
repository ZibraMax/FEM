import matplotlib.pyplot as plt
from FEM.Mesh.Geometry import Geometry

geometria = Geometry.loadGiDMsh('Mesh_tests/GiD/test.msh')
geometria.show()
plt.show()
