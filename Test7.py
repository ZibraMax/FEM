import matplotlib.pyplot as plt
from FEM import Mesh

geometria = Mesh.Geometry.loadGiDMsh('Mesh_tests/GiD/test.msh')
geometria.show()
plt.show()
