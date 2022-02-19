from FEM.Geometry.Geometry import Geometry
from FEM.Geometry.Delaunay import Delaunay
from FEM.Elasticity2D import PlaneStress
from FEM.Utils.polygonal import roundCorner, giveCoordsCircle
import matplotlib.pyplot as plt
import numpy as np

coords = [[0, 0], [4, 0], [4, 0.5], [6, 0.5],
          [6, 2.5], [4, 2.5], [4, 3], [0, 3]]
fillets = [{'start_region': 1, 'end_region': 2, 'r': 0.48, 'n': 5},
           {'start_region': 4, 'end_region': 5, 'r': 0.48, 'n': 5}]
holes = []
radi = 0.5
cent = [2, 1.5]
vert, seg = giveCoordsCircle(cent, radi, n=70)
hole = {'center': cent, 'regions': seg, 'vertices': vert}
holes += [hole]

radi = 0.5/2
cent = [5, 1.5]
vert, seg = giveCoordsCircle(cent, radi, n=40)
hole = {'center': cent, 'regions': seg, 'vertices': vert}
holes += [hole]

params = Delaunay._strdelaunay(constrained=True, delaunay=True, a='0.01', o=2)
geometria = Delaunay(coords, params, nvn=2, fillets=fillets, holes_dict=holes)

cbe = geometria.cbOnHole(0, 0, 1, 3*np.pi/2, np.pi/2)
cbe += geometria.cbOnHole(0, 0, 2, 3*np.pi/2, np.pi/2)
geometria.setCbe(cbe)

# geometria.cbe = geometria.cbFromRegion(7, 0, 1)
# geometria.cbe += geometria.cbFromRegion(7, 0, 2)
# geometria.saveMesh('Mesh_tests/pieza_acero')

E = 29000000
v = 0.26
t = 0.5
p0 = 5000
p = p0/2
# geometria.loadOnRegion(3, lambda s: p)
# for i in range(68, 89):
#     geometria.loadOnRegion(i, lambda s: 10000)

geometria.loadOnHole(1, np.pi/2, 3*np.pi/2, lambda s: p0/(2*np.pi*radi/2))
O = PlaneStress(geometria, E, v, t)
O.geometry.show(draw_bc=True, draw_segs=False)
plt.show()
O.solve()
plt.show()
