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
geometria = Geometry.loadmsh('Mesh_tests/triang_beam.msh')
O = PlaneStress(geometria, E, v, b, fy=lambda x: -gamma)
O.solve()
plt.show()
