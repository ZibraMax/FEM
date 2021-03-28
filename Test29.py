import numpy as np
import matplotlib.pyplot as plt
from FEM.PlaneStrain import PlaneStrain
from FEM.Mesh.Delaunay import Delaunay
from FEM.Utils.polygonal import roundCorner, giveCoordsCircle

E = 21000000.0  # MPa
v = 0.2  # m
gamma = 23.54

b = 200
h = 100
l = 20

vertices = [[0.0, 0.0], [b, 0.0], [b, h], [b/2+2*l, h],
            [b/2+l, h+l], [b/2-l, h+l], [b/2-2*l, h], [0.0, h]]


holes = []
radi = 20
cent = [b/2, h/2]
vert, seg = giveCoordsCircle(cent, radi, n=50)
hole = {'center': cent, 'segments': seg, 'vertices': vert}
holes += [hole]

fillets = [{'start_segment': 2, 'end_segment': 3, 'r': 20, 'n': 10},
           {'start_segment': 3, 'end_segment': 4, 'r': 20, 'n': 10},
           {'start_segment': 4, 'end_segment': 5, 'r': 20, 'n': 10},
           {'start_segment': 5, 'end_segment': 6, 'r': 20, 'n': 10}]

params = Delaunay._strdelaunay(
    constrained=True, delaunay=True, a='7', o=2)
geometria = Delaunay(vertices, params, holes_dict=holes, fillets=fillets)
geometria.show()
plt.show()
O = PlaneStrain(geometria, E, v, fy=lambda x: -gamma)
cb = O.geometry.cbFromSegment(0, 0, 1)
cb += O.geometry.cbFromSegment(0, 0, 2)
cb += O.geometry.cbFromSegment(1, 0, 1)
cb += O.geometry.cbFromSegment(7, 0, 1)
O.cbe = cb
O.solve()
plt.show()
