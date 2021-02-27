import numpy as np
import matplotlib.pyplot as plt
import FEM
from FEM import Mesh

E = 21000000.0 #MPa
v = 0.2 #m
gamma = 23.54
geometria = Mesh.Geometry.loadmsh('Mesh_tests/Talud_fix.msh')
geometria.maskFromSegments()
O = FEM.PlaneStrain(geometria,E,v,fy=lambda x:-gamma*2.5)
O.solve()
plt.show()