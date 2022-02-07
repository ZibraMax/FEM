#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from FEM.Torsion2D import Torsion2D
from FEM.Geometry.Delaunay import Delaunay
from FEM.Geometry.Geometry import Geometry

a = 0.3
b = 0.3
tw = 0.05
tf = 0.05
E = 200000
v = 0.27
G = E / (2 * (1 + v))
phi = 1
vertices = [
    [0, 0],
    [a, 0],
    [a, tf],
    [a / 2 + tw / 2, tf],
    [a / 2 + tw / 2, tf + b],
    [a, tf + b],
    [a, 2 * tf + b],
    [0, 2 * tf + b],
    [0, tf + b],
    [a / 2 - tw / 2, tf + b],
    [a / 2 - tw / 2, tf],
    [0, tf],
]
params = Delaunay._strdelaunay(constrained=True, delaunay=True,
                               a='0.00003', o=2)
geometria = Delaunay(vertices, params)

# geometria.saveMesh('Mesh_tests/I_test')
# geometria = Mesh.Geometry.loadmsh('Mesh_tests/I_test.msh')

print(len(geometria.elements))
O = Torsion2D(geometria, G, phi)
O.solve()
plt.show()
