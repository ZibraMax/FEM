import numpy as np
import matplotlib.pyplot as plt
import FEM
from FEM import Mesh

E = 21000000.0 #MPa
v = 0.2 #m
h = 0.6 #m
b = 0.3 #m
L = 2.5 #m
a = h**2/100
gamma = 23.54
geometria = Mesh.Geometry.loadmsh('Mesh_tests/beam_serendipity.msh')
geometria.cbe = geometria.cbFromSegment(3,0,1)
geometria.cbe += geometria.cbFromSegment(3,0,2)
O = FEM.PlaneStress(geometria,E,v,b,fy=lambda x:-gamma)
O.solve()
plt.show()