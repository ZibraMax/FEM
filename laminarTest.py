from FEM import Geometry2D, PlaneStressSparse, PlaneStressNonLocalSparse, Delaunay, Utils
import numpy as np
import matplotlib.pyplot as plt

from FEM.Geometry.Region import Region1D, Region2D


h = 0.6
b = 0.2
L = 2.5
A = 100


Es = np.array([21.0, 35.0, 28.0, 21.0])*10**6
vs = np.array([0.2, 0.3, 0.25, 0.25])
ncapas = len(Es)
h_capa = h/ncapas

vertices = np.array([
    [0.0, 0.0],
    [L, 0.0],
    [L, h],
    [0.0, h]])


nx = 100
ny = 20
nodes, dicts = Utils.enmalladoFernando(L, h, nx, ny)

bordes = [Region1D(np.array([[0.0, h], [L, h]]))]
for i in range(ncapas):
    capa = Region2D(np.array(
        [[0.0, h_capa*i], [L, h_capa*i], [L, h_capa*(i+1)], [0.0, h_capa*(i+1)]]), desc=f'{Es[i]/10e9:.1f},{vs[i]}')
    bordes += [capa]
bordes += [Region1D(np.array([[0.0, .0], [0.0, h]]))]
bordes += [Region1D(np.array([[L, 0.0], [L, h]]))]

geo = Geometry2D(
    dictionary=dicts,
    gdls=nodes,
    types=['C2V']*len(dicts),
    nvn=2,
    regions=bordes,
    fast=True)
cb = []
cb += geo.generateBCFromCoords(0.0, h/2, 0.0, 1)
cb += geo.generateBCFromCoords(0.0, h/2, 0.0, 2)
cb += geo.generateBCFromCoords(L, h/2, 0.0, 2)

cb += geo.cbFromRegion(-1, 0.0, 2)
cb += geo.cbFromRegion(-2, 0.0, 2)
geo.setCbe(cb)
for i in range(ncapas):
    ele = geo.giveElementsOfRegion(i+1)
    for e in ele:
        e.properties['E'] = Es[i]
        e.properties['v'] = vs[i]
E = np.zeros([len(geo.elements)])
v = np.zeros([len(geo.elements)])

for i in range(len(geo.elements)):
    E[i] = geo.elements[i].properties['E']
    v[i] = geo.elements[i].properties['v']


def carga(s):
    return -A*np.sin(np.pi*s/L)


geo.loadOnRegion(0, fy=carga)
# params = Delaunay._strdelaunay(True, True, L*h/1000, 30, 2)
# geo = Delaunay(vertices=vertices, params=params, nvn=2, fast=True)


geo.show(draw_segs=True, draw_bc=False, label_bc=False)
plt.show()

O = PlaneStressSparse(geo, E, v, b, verbose=True)
O.solve(plot=False)
O.postProcess(levels=30)
O.exportJSON('viga.json')
plt.show()
