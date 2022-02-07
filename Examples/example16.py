import numpy as np
import matplotlib.pyplot as plt
from FEM.PlaneStressNonLocal import PlaneStressNonLocal
from FEM.Geometry.Geometry import Geometry

E = 2.1*10**6
v = 0.2
u = 0.001
a = 5
t = 0.5
l = 0.1
z1 = 0.5
geometria = Geometry.loadmsh('Mesh_tests/EnmalladoTesis.msh')
geometria.generateSegmentsFromCoords([0, 0], [a, 0])
geometria.generateSegmentsFromCoords([a, 0], [a, a])
geometria.generateSegmentsFromCoords([a, a], [0, a])
geometria.generateSegmentsFromCoords([0, a], [0, 0])
cb = geometria.cbFromSegment(3, 0, 1)
cb += geometria.cbFromSegment(3, 0, 2)
cb += geometria.cbFromSegment(1, u, 1)
geometria.setCbe(cb)
geometria.show(draw_bc=True, label_bc=True)
plt.show()


def af(l0, rho):
    return l0*np.exp(-rho)


O = PlaneStressNonLocal(geometria, E, v, t, l, z1, Lr=6*l, af=af)
O.elementMatrices()
O.ensembling()
O.borderConditions()
O.solveES(path='NonLocalGT.csv')
O.postProcess()
plt.show()

_X, U1, U2, U3, U = O.profile([0, 0.019], [a, 0.019])
np.savetxt('a2gt.csv', [_X, U1], delimiter=',')
plt.show()
