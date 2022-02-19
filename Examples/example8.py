import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity2D import PlaneStress
from FEM.Geometry.Geometry import Geometry

a = 5
u0 = 0.001

E = 2.1*10**6
v = 0.2
t = 0.5

geometria = Geometry.loadmsh("Mesh_tests/rect.msh")
O = PlaneStress(geometria, E, v, t)
cbe = O.geometry.cbFromRegion(3, 0, 1)
cbe += O.geometry.cbFromRegion(3, 0, 2)
cbe += O.geometry.cbFromRegion(1, u0, 1)
O.cbe = cbe
O.solve()
plt.show()
