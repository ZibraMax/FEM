from FEM.Heat2D import Heat2D
from FEM.Mesh.Delaunay import Delaunay
from FEM.Utils.polygonal import roundCorner, giveCoordsCircle
import numpy as np
import matplotlib.pyplot as plt
# https://github.com/Samson-Mano/2D_Heat_transfer
k = 2
t0 = 140

beta = 1.5
Ta = 20

radi = 4
cent = [0, 0]
vert_orig, seg_orig = giveCoordsCircle(cent, radi, n=150)

vert, seg = giveCoordsCircle([0, 1], 1, n=50)
holes = [{'center': [0, 1], 'segments': seg, 'vertices': vert}]
params = Delaunay._strdelaunay(constrained=True, delaunay=True, a='0.03', o=2)
geometria = Delaunay(vert_orig, params, 1, holes_dict=holes)
geometria.setCbe(geometria.cbOnHole(0, t0, 1))
O = Heat2D(geometria, k, k)
for i in range(len(seg_orig)):
    O.defineConvectiveBoderConditions(i, beta, Ta)
O.geometry.show()
plt.show()
O.solve()
plt.show()
