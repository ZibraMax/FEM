import numpy as np
import matplotlib.pyplot as plt
from FEM.PlaneStressNonLocal import PlaneStressNonLocal
from FEM.Mesh.Geometry import Geometry

E = 21000000.0  # MPa
v = 0.2  # m
h = 0.6  # m
b = 0.3  # m
L = 2.5  # m
a = h**2/100
gamma = 23.54

z1 = 0.17
l = 0.05
Lr = 6*l

q0 = -gamma*b
EI = E*b*h**3/12
k = l


def af(l0, rho):
    return 2*np.pi/Lr/Lr*(1-(rho)**2/(Lr)**2)*(rho <= Lr)


def w(x):
    return q0/EI*(x**4/24-L*x**3/6+L**2*x**2/4)+q0*k/2/EI*(-k*x**2+2*k*L*x+L**2*x)


geometria = Geometry.loadmsh('Examples/No Local/viga_empotrada.msh')
geometria.generateSegmentsFromCoords([0, 0], [L, 0])
geometria.generateSegmentsFromCoords([L, 0], [L, h])
geometria.generateSegmentsFromCoords([L, h], [0, h])
geometria.generateSegmentsFromCoords([0, h], [0, 0])
cbe = geometria.cbFromSegment(3, 0, 1)
cbe += geometria.cbFromSegment(3, 0, 2)
geometria.cbe = cbe
O = PlaneStressNonLocal(geometria, E, v, b, l, z1,
                        Lr, af, fy=lambda x: -gamma*b)
O.geometry.mask = None
O.solve()
# O.solveFromFile('U_1.csv', plot=False)


_X, U1, U2, U3, UU = O.profile([0, h*0.49], [L, h*0.5], n=1000)
# plt.savefig('local.svg', transparent=True)
plt.show()

# plt.plot(_X, U1, '-', c='k', label='No Local')
# plt.plot(_X2, U12, '--', c='k', label='Local')
# plt.grid()
# plt.xlabel('x')
# plt.ylabel(r'$\varepsilon_x$')
# plt.legend()
# plt.show()

# plt.plot(_X, U2, '-', c='k', label='No Local')
# plt.plot(_X2, U22, '--', c='k', label='Local')
# plt.grid()
# plt.xlabel('x')
# plt.ylabel(r'$\varepsilon_y$')
# plt.legend()
# plt.show()

# plt.plot(_X, U3, '-', c='k', label='No Local')
# plt.plot(_X2, U32, '--', c='k', label='Local')
# plt.grid()
# plt.xlabel('x')
# plt.ylabel(r'$\varepsilon_{xy}$')
# plt.legend()
# plt.show()

plt.plot(_X, np.array(UU)[:, 1, 0], '-', c='k', label='No Local')
plt.plot(_X, w(np.array(_X)), '-.', c='k', label='No Local Analítica')
plt.grid()
plt.xlabel('x')
plt.ylabel(r'$\Delta$')
plt.legend()
plt.show()

# plt.savefig('nolocal.svg', transparent=True)
