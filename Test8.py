import numpy as np
import matplotlib.pyplot as plt
import FEM
from FEM import Mesh

a = 5
u0 = 0.001

E = 2.1*10**6
v = 0.2
t = 0.5

geometria = Mesh.Geometry.loadmsh("Mesh_tests/Square_torsion.msh")
O = FEM.PlaneStress(geometria,E,v,t)
cbe = O.geometry.cbFromSegment(3,0,1)
cbe += O.geometry.cbFromSegment(3,0,2)
cbe += O.geometry.cbFromSegment(1,u0,1)
O.cbe = cbe
O.solve()
plt.show()