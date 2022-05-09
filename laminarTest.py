from FEM import Geometry2D, PlaneStressSparse, PlaneStressNonLocalSparse, Delaunay, Utils
import numpy as np
import matplotlib.pyplot as plt

from FEM.Geometry.Region import Region1D

ncapas = 4

h = 0.6
b = 0.2
L = 2.5
A = 1000

h_capa = h/ncapas

Es = np.array([21, 35, 28])*10**9
vs = np.array([0.2, 0.3, 0.25])


vertices = np.array([
    [0.0, 0.0],
    [L, 0.0],
    [L, h],
    [0.0, h]])


nx = 100
ny = 20
nodes, dicts = Utils.enmalladoFernando(L, h, nx, ny)

bordes = [Region1D(np.array([[0.0, h], [L, h]]))]

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
geo.setCbe(cb)


def carga(s):
    return A*np.sin(np.pi*s/L)


geo.loadOnRegion(0, fy=carga)
# params = Delaunay._strdelaunay(True, True, L*h/1000, 30, 2)
# geo = Delaunay(vertices=vertices, params=params, nvn=2, fast=True)


# geo.show(draw_bc=True, label_bc=True)
# plt.show()

O = PlaneStressSparse(geo, 21*10**9, 0.2, b, verbose=True)
O.solve(mult=10000)
O.exportJSON('viga.json')
plt.show()
