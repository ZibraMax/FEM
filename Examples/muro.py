import numpy as np
import matplotlib.pyplot as plt
from FEM import PlaneStressSparse
from FEM import Delaunay
from FEM.Geometry.Region import Region1D

b = 5.65
h = 5.7
u0 = -0.03
vertices = [[0.0, 0.0], [b, 0.0], [b, h], [0.0, h]]
params = Delaunay._strdelaunay(True, True, 0.01, 30, 2)

b_p = 5.24-4.03
h_p = 2.53
c_puerta = np.array([[0.0, 0.0], [b_p, 0.0], [b_p, h_p], [0.0, h_p]])
c_puerta[:, 0] += b/2+1/2*b_p

holes = []
puerta = {'center': np.average(
    c_puerta, axis=0), 'regions': [[0, 1], [1, 2], [2, 3], [3, 0]], 'vertices': c_puerta}
holes += [puerta]

c_ventana = np.array([[0.0, 0.0], [b_p, 0.0], [b_p, b_p], [0.0, b_p]])
c_ventana[:, 0] += b/2-3/2*b_p
c_ventana[:, 1] += h-(5.81-5.24)-b_p

ventana = {'center': np.average(
    c_ventana, axis=0), 'regions': [[0, 1], [1, 2], [2, 3], [3, 0]], 'vertices': c_ventana}
holes += [ventana]

c_ventana = np.array([[0.0, 0.0], [b_p, 0.0], [b_p, b_p], [0.0, b_p]])
c_ventana[:, 0] += b/2+1/2*b_p
c_ventana[:, 1] += h-(5.81-5.24)-b_p

ventana = {'center': np.average(
    c_ventana, axis=0), 'regions': [[0, 1], [1, 2], [2, 3], [3, 0]], 'vertices': c_ventana}
holes += [ventana]


geo = Delaunay(vertices=vertices, params=params,
               nvn=2, holes_dict=holes, fast=True)
h_apoyo = 5.7-(2.65+2.69)
reg = Region1D(np.array([[b, 2.65-h_apoyo], [b, 2.56+h_apoyo]]))
geo.addRegions([reg])
cb = geo.cbFromRegion(16, u0, 1)
cb += geo.cbFromRegion(0, 0.0, 1)
cb += geo.cbFromRegion(0, 0.0, 2)
geo.setCbe(cb)

# geo.show()
# plt.show()

E = 97000.0
v = 0.3
t = 0.65

O = PlaneStressSparse(geo, E, v, t, verbose=True)
O.solve(mult=100/3)
O.exportJSON("muro.json")
plt.show()
