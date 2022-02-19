from FEM.Geometry import Delaunay
from FEM.Torsion2D import Torsion2D
import matplotlib.pyplot as plt

coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
fillets = [{'start_region': 3, 'end_region': 0, 'r': 0.5, 'n': 10}]
params = Delaunay._strdelaunay(constrained=True, delaunay=True, a='0.005', o=2)
geometria = Delaunay(coords, params, nvn=1, fillets=fillets)
geometria.show()
plt.show()
O = Torsion2D(geometria, 1, 1)
O.geometry.mask = None
O.solve()
plt.show()
