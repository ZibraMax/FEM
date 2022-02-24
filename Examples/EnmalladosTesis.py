import numpy as np
import tetgen
import pyvista as pv
from FEM.Geometry import Geometry3D
pv.set_plot_theme('document')
# sphere = pv.Sphere()
# tet = tetgen.TetGen(sphere)
# coords, dicc = tet.tetrahedralize(order=1, mindihedral=20, minratio=1.5)
# grid = tet.grid
# grid.plot(show_edges=True)
# geo = Geometry3D(dicc.tolist(), coords, ['TE1V']*len(dicc), nvn=3, fast=True)
# geo.exportJSON('SPHERE.json')

ledge = 10.0
height = ledge*6**0.5/3

v = np.array([[0, 0, 0], [ledge, 0, 0.0],
              [0, ledge, 0], [0, 0, height]])
f = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
tgen = tetgen.TetGen(v, f)
coords, dicc = tgen.tetrahedralize(switches='-pqa10i')
tgen.plot(show_edges=True)
