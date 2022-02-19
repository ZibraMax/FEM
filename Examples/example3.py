#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from FEM.Torsion2D import Torsion2D
from FEM.Geometry import Geometry2D, Delaunay

a = 1
b = 1
E = 200000
v = 0.27
G = E/(2*(1+v))*0+1
phi = 1

# coords = np.array([[0, 0], [1, 0], [1, 1], [0, 1.0]])
# params = Delaunay._strdelaunay(a=0.0001, o=2)
# geo = Delaunay(coords, params)
# geo.exportJSON('Examples/Mesh_tests/Square_torsion.json')
geometria = Geometry2D.importJSON('Examples/Mesh_tests/Square_torsion.json')
# geometria.show()
# plt.show()
O = Torsion2D(geometria, G, phi, verbose=True)
O.solve()
plt.show()
