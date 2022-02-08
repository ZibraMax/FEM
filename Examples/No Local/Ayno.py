import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity2D import PlaneStressNonLocal
from FEM.Geometry.Delaunay import Delaunay

E = 2.1*10**6
v = 0.2
u = 0.001
t = 0.5
l = 0.1
z1 = 0.5
LR = 6*l
fact = 10
P = 1
a = 5
coords = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])*a
params = Delaunay._strdelaunay(constrained=True, delaunay=True, a='0.02', o=2)
geometria = Delaunay(coords, params, nvn=2)
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
O.solveES()

O.postProcess()
plt.show()

_X, U1, U2, U3 = O.profile([0, 0.019], [a, 0.019])
np.savetxt('a2.csv', [_X, U1], delimiter=',')
plt.show()
