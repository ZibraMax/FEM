import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity2D import PlaneStress
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
geometria.generateRegionFromCoords([0, 0], [L, 0])
geometria.generateRegionFromCoords([L, 0], [L, h])
geometria.generateRegionFromCoords([L, h], [0, h])
geometria.generateRegionFromCoords([0, h], [0, 0])
cbe = geometria.cbFromRegion(3, 0, 1)
cbe += geometria.cbFromRegion(3, 0, 2)
cbe += geometria.cbFromRegion(1, 0, 1)
cbe += geometria.cbFromRegion(1, 0, 2)
geometria.cbe = cbe
O = PlaneStress(geometria, E, v, b, fy=lambda x: -gamma*b)

O.geometry.mask = None
O.solve()
np.savetxt('U_2_Local.csv', O.U, delimiter=',')
plt.show()
