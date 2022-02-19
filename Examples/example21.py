import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity2D import PlaneStress
from FEM.Geometry.Geometry import Geometry

b = 120
h = 160
t = 0.036
E = 30*10**(6)
v = 0.25

gdls = [[0, 0], [b, 0], [0, h], [b, h]]
elemento1 = [0, 1, 3, 2]
dicc = [elemento1]
tipos = ['C1V']
segmentos = [[0, 1], [1, 3], [3, 2], [2, 0]]
geometria = Geometry(dicc, gdls, tipos, nvn=2, regions=segmentos)
cbe = geometria.cbFromRegion(3, 0, 1)
cbe += geometria.cbFromRegion(3, 0, 2)
geometria.cbe = cbe
geometria.loadOnRegion(1, fx=lambda s: 10, fy=lambda s: 0)

O = PlaneStress(geometria, E, v, t)
O.elementMatrices()
O.ensembling()
O.borderConditions()
O.solveES()
O.postProcess()
# plt.show()

print(O.giveStressPoint(np.array([[60], [80]])))
