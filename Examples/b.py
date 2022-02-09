import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity3D import Elasticity
from FEM.Geometry.Geometry import Geometry


E = 21000000.0
v = 0.2
h = 0.6
b = 0.3
L = 2.5
gamma = 23.54

_a = L
_b = h
_c = b

nx = 50
ny = 8
nz = 7

dx = _a/nx
dy = _b/ny
dz = _c/nz

coords = []

for i in range(nx+1):
    x = i*dx
    for j in range(ny+1):
        y = j*dy
        for k in range(nz+1):
            z = k*dz
            coords += [[x, y, z]]

dicc = []


def node(i, j, k): return i*(ny+1)*(nz+1)+j*(nz+1)+k


for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            node1 = node(i, j, k)
            node2 = node(i+1, j, k)
            node3 = node(i+1, j+1, k)
            node4 = node(i, j+1, k)

            node5 = node(i, j, k+1)
            node6 = node(i+1, j, k+1)
            node7 = node(i+1, j+1, k+1)
            node8 = node(i, j+1, k+1)

            dicc += [[node1, node2, node3, node4, node5, node6, node7, node8]]


def fy(x): return -gamma


geometria = Geometry(dicc, coords, ["B1V"]*len(dicc), nvn=3, fast=True)
cbe = []
for i in range(len(coords)):
    if 0.0 == coords[i][0]:
        cbe += [[i*3, 0.0]]
        cbe += [[i*3+1, 0.0]]
        cbe += [[i*3+2, 0.0]]
geometria.cbe = cbe

O = Elasticity(geometria, E, v, gamma, fy=fy, verbose=True)
O.solve()
coords = [[0.0, 0.0, b/2], [L, 0.0, b/2], [L, h, b/2], [0, h, b/2]]
x, resultU, resultDU = O.profile(coords, 20)
es = []
e = np.zeros([3, 3])
for res in resultDU:
    e = []
    defo = res[0]
    ex = defo[0, 0]
    ey = defo[1, 1]
    ez = defo[2, 2]
    exy = (defo[0, 1]+defo[1, 0])
    exz = (defo[0, 2]+defo[2, 0])
    eyz = (defo[1, 2]+defo[2, 1])
    e.append(ex)
    e.append(ey)
    e.append(ez)
    e.append(exy)
    e.append(exz)
    e.append(eyz)
    es += [e]

C = E/((1.0+v)*(1.0-2.0*v))*np.array([
    [1.0-v, v, v, 0.0, 0.0, 0.0],
    [v, 1.0-v, v, 0.0, 0.0, 0.0],
    [v, v, 1.0-v, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, (1.0-2.0*v)/2.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, (1.0-2.0*v)/2.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, (1.0-2.0*v)/2.0]])


es = np.array(es)
sig = es@C
X = x[:, 0]
Y = x[:, 1]
Z = sig[:, 0]
surf = plt.tricontourf(X, Y, Z, cmap='rainbow')
plt.colorbar(surf)
plt.show()
a = 0
