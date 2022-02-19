from operator import ge
import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity2D import PlaneStressNonLocalSparse
from FEM.Geometry.Geometry import Geometry


E = 2.1*10**6
v = 0.2
u = 0.001
t = 0.5
l = 0.1
z1 = 0.5
LR = 6*l
P = 1
a = 5


nx = 30
ny = 30
nz = 0

_a = a
_b = a
_c = t

dx = _a/nx
dy = _b/ny
dz = 0.0

coords = []

for i in range(nx+1):
    x = i*dx
    for j in range(ny+1):
        y = j*dy
        coords += [[x, y]]

dicc = []


def node(i, j, k): return i*(ny+1)*(nz+1)+j*(nz+1)+k
def node3d(i, j, k): return i*(30+1)*(6+1)+j*(6+1)+k


gdls = []

for i in range(nx):
    for j in range(ny):
        k = 0
        node1 = node(i, j, k)
        node2 = node(i+1, j, k)
        node3 = node(i+1, j+1, k)
        node4 = node(i, j+1, k)

        dicc += [[node1, node2, node3, node4]]
for i in range(nx+1):
    for j in range(ny+1):
        k = 3
        gdls += [node3d(i, j, k)*3]
        gdls += [node3d(i, j, k)*3+1]


geometria2D = Geometry(dicc, coords, ['C1V']*len(dicc), 2, fast=True)
# geometria2D.show(draw_labels=True)
# plt.show()
# geometria3D = Geometry.importJSON('3DNonLocal.json', fast=True)

# coords = []
# elements = []
# coords = np.array(geometria3D.gdls)
# coords = coords[coords[:, 2] == 0.25]
# geometria = Geometry.loadmsh(
#     'Examples/Mesh_tests/EnmalladoTesis.msh', fast=True)
# geometria.generateRegionFromCoords([0, 0], [a, 0])
# geometria.generateRegionFromCoords([a, 0], [a, a])
# geometria.generateRegionFromCoords([a, a], [0, a])
# geometria.generateRegionFromCoords([0, a], [0, 0])
# cb = geometria.cbFromSegment(3, 0, 1)
# cb += geometria.cbFromSegment(3, 0, 2)
# cb += geometria.cbFromSegment(1, u, 1)
# geometria.setCbe(cb)


# geometria.show(draw_bc=True, label_bc=True)
# plt.show()


def af(l0, rho):
    return l0*np.exp(-rho)


O = PlaneStressNonLocalSparse(
    geometria2D, E, v, t, l, z1, Lr=6*l, af=af, verbose=True)

arr = np.loadtxt('a_gt.csv')

O.U = arr[np.ix_(gdls)]
for e in O.elements:
    e.setUe(O.U)
O.postProcess()
plt.show()

_X, U1, U2, U3, U = O.profile([0, 0.019], [a, 0.019])
np.savetxt('a2.csv', [_X, U1], delimiter=',')
plt.show()
