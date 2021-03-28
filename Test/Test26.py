from FEM.Mesh.Geometry import Geometry
from FEM.Mesh.Delaunay import Delaunay
from FEM.Torsion2D import Torsion2D
from FEM.Utils.polygonal import roundCorner, giveCoordsCircle
import matplotlib.pyplot as plt
import numpy as np

coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
fillets = [{'start_segment': 3, 'end_segment': 0, 'r': 0.5, 'n': 10}]
params = Delaunay._strdelaunay(constrained=True, delaunay=True, a='0.005', o=2)
geometria = Delaunay(coords, params, nvn=1, fillets=fillets)
geometria.mask = None
geometria.show()
plt.show()
O = Torsion2D(geometria, 1, 1)
O.solve()
plt.show()
