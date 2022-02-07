import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity2D import PlaneStress
from FEM.Geometry.Geometry import Geometry
from FEM.Geometry.Delaunay import Delaunay

b = 120
h = 160
t = 0.036
E = 30*10**(6)
v = 0.25

gdls = [[0, 0], [b, 0],  [b, h], [0, h]]
params = Delaunay._strdelaunay(constrained=True, delaunay=True, a='10', o=2)
geometria = Delaunay(gdls, params, nvn=2)
cbe = geometria.cbFromSegment(3, 0, 1)
cbe += geometria.cbFromSegment(3, 0, 2)
geometria.cbe = cbe
geometria.loadOnSegment(1, fx=lambda s: 10, fy=lambda s: 0)
O = PlaneStress(geometria, E, v, t)
O.elementMatrices()
O.ensembling()
O.borderConditions()
O.solveES()
O.postProcess()
plt.show()

print(O.giveStressPoint(np.array([[60], [80]])))
