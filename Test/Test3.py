#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import FEM
from FEM.Mesh import *

a = 1
b = 1
E = 200000
v = 0.27
G = E/(2*(1+v))*0+1
phi = 1
geometria = Geometry.loadmsh('Mesh_tests/Square_torsion.msh')
# geometria.show()
# plt.show()
O = FEM.Torsion2D(geometria, G, phi)
O.solve()
plt.show()
