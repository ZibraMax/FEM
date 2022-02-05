import scipy
import json
import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity3D import Elasticity
from FEM.Mesh.Geometry import Geometry


E = 21000000.0
v = 0.2
h = 1
b = 1
L = 1

gamma = 23.54

_a = L
_b = h
_c = b

nx = 10
ny = 10
nz = 10

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


geometria = Geometry(dicc, coords, ["B1V"]*len(dicc), nvn=3, fast=True)

O = Elasticity(geometria, E, v, gamma, verbose=True)
O.solve()
omh = scipy.sparse.linalg.spsolve(O.M, O.K)
eiv, eigvec = scipy.sparse.linalg.eigs(omh, k=1, which='LM')
y = O.geometry.exportJSON()
pjson = json.loads(y)
pjson["disp_field"] = eigvec.real.flatten().tolist()
y = json.dumps(pjson)
with open("../FEM C++/docs/resources/exported.json", "w") as f:
    f.write(y)
np.savetxt('a.csv', O.U, delimiter=',', fmt='%s')
a = 0
