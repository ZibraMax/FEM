import numpy as np
import matplotlib.pyplot as plt
from FEM.Elasticity2D import PlaneStress
from FEM.Mesh.Geometry import Geometry
from FEM.Utils import enmalladoFernando

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
filename = 'viga_empotrada.msh'
enmalladoFernando(L, h, 80, 20, filename)


def af(l0, rho):
    return 2*np.pi/Lr/Lr*(1-(rho)**2/(Lr)**2)*(rho <= Lr)


geometria = Geometry.loadmsh(filename)
geometria.generateSegmentsFromCoords([0, 0], [L, 0])
geometria.generateSegmentsFromCoords([L, 0], [L, h])
geometria.generateSegmentsFromCoords([L, h], [0, h])
geometria.generateSegmentsFromCoords([0, h], [0, 0])
cbe = geometria.cbFromSegment(3, 0, 1)
cbe += geometria.cbFromSegment(3, 0, 2)
geometria.cbe = cbe
O = PlaneStress(geometria, E, v, b, fy=lambda x: -gamma*b)
O.geometry.mask = None
O.solve()
np.savetxt('U_1_Local.csv', O.U, delimiter=',')
plt.show()
