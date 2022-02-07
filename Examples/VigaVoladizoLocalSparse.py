import scipy
from fileinput import filename
import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity2D import PlaneStressSparse
from FEM.Geometry.Geometry import Geometry
from FEM.Utils.polygonal import enmalladoFernando

E = 30000.0*6895.0  # KPa
v = 0.25
h = 2*2.54/100  # m
b = 1*2.54/100  # m
L = 10*2.54/100  # m
# a = h**2/100
t0 = 1034.21  # kpa

rho = 7.860

filename = "viga_voladizo.msh"
nx = 500
ny = 10
enmalladoFernando(L, h, nx, ny, filename)
geometria = Geometry.loadmsh(filename, fast=True)
geometria.generateSegmentsFromCoords([0, 0], [L, 0])
geometria.generateSegmentsFromCoords([L, 0], [L, h])
geometria.generateSegmentsFromCoords([L, h], [0, h])
geometria.generateSegmentsFromCoords([0, h], [0, 0])

cbe = geometria.cbFromSegment(1, 0.0, 1)
cbe += geometria.generateBCFromCoords(L, h/2, 0.0, 2)

geometria.loadOnSegment(3, fy=lambda s: -t0*b)

geometria.setCbe(cbe)
O = PlaneStressSparse(geometria, E, v, b, rho=rho, verbose=True)
O.geometry.mask = None
O.solve()
# np.savetxt('U_2_Local.csv', O.U, delimiter=',')
plt.show()

# O.elementMatrices()
# O.ensembling()
# O.borderConditions()
# OM = scipy.sparse.linalg.spsolve(O.M, O.K)
# omega = scipy.sparse.linalg.eigs(
#     OM, 10, return_eigenvectors=False, which="SM")**0.5
# print(omega)
