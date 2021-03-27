from FEM.Mesh.Geometry import Geometry
from FEM.Mesh.Delaunay import Delaunay
from FEM.PlaneStrain import PlaneStrain
from FEM.Utils.polygonal import roundCorner, giveCoordsCircle
import matplotlib.pyplot as plt
import numpy as np

coords = [[0, 0], [1, 0], [1, 1], [0, 1]]
fillets = [{'start_segment': 0, 'end_segment': 1, 'r': 0.1, 'n': 5}]
params = Delaunay._strdelaunay(constrained=True, delaunay=True, a='0.01', o=2)
geometria = Delaunay(coords, params, nvn=1, fillets=fillets)
geometria.show()
plt.show()
