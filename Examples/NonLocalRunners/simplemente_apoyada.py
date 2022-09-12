from FEM import Geometry2D, PlaneStressNonLocalSparse, PlaneStressSparse, Utils
from FEM.Elasticity2D import PlaneStress
from FEM.Solvers import LinealEigen
import numpy as np
import matplotlib.pyplot as plt
from FEM.Geometry.Region import Region1D, Region2D


h = 0.6
b = 0.2
L = 10
E = 30e6
rho = 1
v = 0.2
l = 0.5
z1 = 1.0
Lr = 6*l


vertices = np.array([
    [0.0, 0.0],
    [L, 0.0],
    [L, h],
    [0.0, h]])


nx = 100
ny = 8
nodes, dicts = Utils.enmalladoFernando(L, h, nx, ny)

bordes = [Region1D(np.array([[0.0, h], [L, h]]))]  # Superior
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
# cb += geo.generateBCFromCoords(0.0, h/2, 0.0, 1)
# cb += geo.generateBCFromCoords(L, h/2, 0.0, 1)

cb += geo.cbFromRegion(1, 0.0, 1)
cb += geo.cbFromRegion(1, 0.0, 2)
geo.setCbe(cb)


def af(l0, rho):
    return l0*np.exp(-rho)


# def carga(s):
#     return -A*np.sin(np.pi*s/L)


# geo.loadOnRegion(0, fy=carga)
# params = Delauqnay._strdelaunay(True, True, L*h/1000, 30, 2)
# geo = Delaunay(vertices=vertices, params=params, nvn=2, fast=True)


geo.show(draw_segs=True, draw_bc=True, label_bc=False)
plt.show()

# O = PlaneStressNonLocalSparse(
#     geo, E, v, b, l, z1, Lr, af, rho, solver=LinealEigen, verbose=True)
O = PlaneStressSparse(
    geo, E, v, b, rho, solver=LinealEigen, verbose=True)
O.solve(plot=False)

O.exportJSON('viga.json')
plt.show()
