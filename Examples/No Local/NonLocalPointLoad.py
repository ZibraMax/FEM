import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity2D import PlaneStressNonLocal
from FEM.Geometry.Delaunay import Delaunay

E = 2.1*10**6
v = 0.2
u = 0.001
t = 1
l = 0.1
z1 = 0.5
LR = 6*l
fact = 10
P = 1
coords = np.array([[0, 0], [2*LR, 0], [2*LR, 10*LR], [0, 10*LR]])
params = Delaunay._strdelaunay(constrained=True, delaunay=True, a='0.003', o=1)
geometria = Delaunay(coords, params, nvn=2)
cb = geometria.cbFromSegment(3, 0, 1)
cb += geometria.cbFromSegment(1, 0, 1)
cb += geometria.cbFromSegment(0, 0, 1)
cb += geometria.cbFromSegment(0, 0, 2)
geometria.setCbe(cb)

cb = geometria.generateBCFromCoords(LR, 10*LR, P, nv=2)
geometria.cbn = cb

geometria.show(draw_bc=True, label_bc=True)
plt.show()


def af(l0, rho):
    return 2*np.pi/LR/LR*(1-(rho)**2/(LR)**2)*(rho <= LR)


O = PlaneStressNonLocal(geometria, E, v, t, l, z1, Lr=6*l, af=af)
O.elementMatrices()
O.ensembling()
O.borderConditions()
O.solveES()

O.postProcess()
plt.show()

# _X, U1, U2, U3 = O.profile([0, 0.019], [a, 0.019])
X, Y = O.profileStress([LR, 10*LR], [LR, 0], n=300)
plt.show()
