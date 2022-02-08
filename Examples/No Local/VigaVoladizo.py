import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity2D import PlaneStressNonLocal
from FEM.Geometry.Geometry import Geometry

E = 21000000.0  # MPa
v = 0.2  # m
h = 0.6  # m
b = 0.3  # m
L = 2.5  # m
a = h**2/100
gamma = 23.54

z1 = 0.5
l = 0.05
Lr = 6*l


def af(l0, rho):
    return 2*np.pi/Lr/Lr*(1-(rho)**2/(Lr)**2)*(rho <= Lr)


geometria = Geometry.loadmsh('viga_empotrada.msh')
geometria.generateSegmentsFromCoords([0, 0], [L, 0])
geometria.generateSegmentsFromCoords([L, 0], [L, h])
geometria.generateSegmentsFromCoords([L, h], [0, h])
geometria.generateSegmentsFromCoords([0, h], [0, 0])
cbe = geometria.cbFromSegment(3, 0, 1)
cbe += geometria.cbFromSegment(3, 0, 2)
cbe += geometria.cbFromSegment(1, 0, 1)
cbe += geometria.cbFromSegment(1, 0, 2)
geometria.cbe = cbe
O = PlaneStressNonLocal(geometria, E, v, b, l, z1,
                        Lr, af, fy=lambda x: -gamma*b)
O.geometry.mask = None
O.solve()
np.savetxt('U_2.csv', O.U, delimiter=',')
plt.show()
