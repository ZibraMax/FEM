import json
import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity3D import NonLocalElasticity
from FEM.Geometry.Geometry import Geometry


E = 2.1*10**6
v = 0.2
t = 0.5
l = 0.1
z1 = 0.5
Lr = 6*l
a = 5.0

u0 = 0.001

gamma = 23.54/9.81


def af(rho):
    return (1/(8*np.pi*l**3))*np.exp(-rho)  # No referencia, sacada a mano


nx = 30
ny = 30
nz = 6

_a = a
_b = a
_c = t

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
cbe = []
for i in range(len(coords)):
    if 0.0 == coords[i][0]:
        cbe += [[i*3, 0.0]]
        cbe += [[i*3+1, 0.0]]
        cbe += [[i*3+2, 0.0]]

for i in range(len(coords)):  # ATENCIÓN SE ESTAN COMPARANDO FLOTANTES TENER MUCHO CUIDADO!!!!
    if a == coords[i][0]:
        cbe += [[i*3, u0]]

geometria.cbe = cbe
geometria.exportJSON('3DNonLocal.json')
O = NonLocalElasticity(geometria, E, v, gamma, l, z1, Lr, af, verbose=True)
arr = []
for e in O.elements:
    arr += [e.enl]
# np.savetxt('nolocales.csv', arr, delimiter=',')
O.solve()
y = O.geometry.exportJSON()
pjson = json.loads(y)
pjson["disp_field"] = O.U.tolist()
y = json.dumps(pjson)
with open("exported.json", "w") as f:
    f.write(y)
np.savetxt('a.csv', O.U, delimiter=',', fmt='%s')
a = 0
