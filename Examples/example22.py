import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity2D import PlaneStressSparse
from FEM.Geometry import Delaunay

b = 120
h = 160
t = 0.036
E = 30*10**(6)
v = 0.25

gdls = [[0, 0], [b, 0],  [b, h], [0, h]]
params = Delaunay._strdelaunay(constrained=True, delaunay=True, a='10', o=2)
geometria = Delaunay(gdls, params, nvn=2, fast=True)
cbe = geometria.cbFromRegion(3, 0, 1)
cbe += geometria.cbFromRegion(3, 0, 2)
geometria.cbe = cbe
geometria.loadOnRegion(1, fx=lambda s: 10, fy=lambda s: 0)
O = PlaneStressSparse(geometria, E, v, t, verbose=True)
O.elementMatrices()
O.ensembling()
O.borderConditions()
O.solveES()
O.postProcess()
plt.show()

print(O.giveStressPoint(np.array([[60], [80]])))
