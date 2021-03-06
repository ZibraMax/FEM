import numpy as np
import matplotlib.pyplot as plt
from FEM.Mesh.Geometry import Geometry
from FEM.PlaneStrain import PlaneStrain

E = 21000000.0  # MPa
v = 0.2  # m
gamma = 23.54
geometria = Geometry.loadmsh('Mesh_tests/Talud_fix.msh')
geometria.maskFromSegments()
print(geometria.mask)
O = PlaneStrain(geometria, E, v, fy=lambda x: -gamma*2.5)
O.solve()
plt.show()
