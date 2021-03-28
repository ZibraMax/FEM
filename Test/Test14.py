import numpy as np
import matplotlib.pyplot as plt
from FEM.PlaneStress import PlaneStress
from FEM.Mesh.Geometry import Geometry

E = 21000000.0  # MPa
v = 0.2  # m
h = 0.6  # m
b = 0.3  # m
L = 2.5  # m
a = h**2/100
gamma = 23.54
geometria = Geometry.loadmsh('Mesh_tests/beam_serendipity.msh')
geometria.cbe = geometria.cbFromSegment(3, 0, 1)
geometria.cbe += geometria.cbFromSegment(3, 0, 2)
O = PlaneStress(geometria, E, v, b, fy=lambda x: -gamma)
O.solve('Fast_test_ser.csv')
plt.show()
