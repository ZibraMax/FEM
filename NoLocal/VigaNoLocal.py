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

z1 = 0.5
l = 0.05
Lr = 6*l

q0 = -gamma*b
EI = E*b*h**3/12
tao = l


def af(l0, rho):
    return 2*np.pi/Lr/Lr*(1-(rho)**2/(Lr)**2)*(rho <= Lr)


sinh = np.sinh
cosh = np.cosh
sqrtz = z1**0.5


def w3(x):
    return q0*x**4/24/EI-q0*L*x**3/12/EI+q0*x**2*((L**3+12*L*(z1-1)*tao**2+24*(z1-1)*tao**3)*sinh(L/2/sqrtz/tao)+L**3*sqrtz*cosh(L/2/sqrtz/tao))/24/EI/((L-2*(z1-1)*tao)*sinh(L/2/sqrtz/tao)+L*sqrtz*cosh(L/2/sqrtz/tao))-(q0*L*(z1-1)*tao*(L**2+6*L*tao+12*tao**2)*sinh(L/2/sqrtz/tao)*x)/12/EI/((L-2*(z1-1)*tao)*sinh(L/2/sqrtz/tao)+L*sqrtz*cosh(L/2/sqrtz/tao))+(q0*L*(z1-1)*sqrtz*tao**2*(L**2+6*L*tao+12*tao**2)*sinh(x/2/sqrtz/tao)*sinh((L-x)/(2*sqrtz*tao)))/6/EI/((L-2*(z1-1)*tao)*sinh(L/2/sqrtz/tao)+L*sqrtz*cosh(L/2/sqrtz/tao))


def w2(x):
    return q0*x**4/24/EI-q0*L*x**3/12/EI+q0*x**2*((L**3+12*L*(z1-1)*tao**2+24*(z1-1)*tao**3))/24/EI/((L+2*(1-sqrtz)*tao))+(q0*L*(1-sqrtz)*(L**2+6*L*tao+12*tao**2)*tao*x)/12/EI/((L+2*(1-sqrtz)*tao))-(q0*L*sqrtz*(1-sqrtz)*(L**2+6*L*tao+12*tao**2)*tao**2)/12/EI/((L+2*(1-sqrtz)*tao))*(1-np.exp(-x/sqrtz/tao)-np.exp(-(L-x)/sqrtz/tao))


geometria = Geometry.loadmsh('Examples/No Local/viga_empotrada.msh')
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
O.solveFromFile('U_2.csv', plot=False)


geometria = Geometry.loadmsh('Examples/No Local/viga_empotrada.msh')
geometria.generateSegmentsFromCoords([0, 0], [L, 0])
geometria.generateSegmentsFromCoords([L, 0], [L, h])
geometria.generateSegmentsFromCoords([L, h], [0, h])
geometria.generateSegmentsFromCoords([0, h], [0, 0])
cbe = geometria.cbFromSegment(3, 0, 1)
cbe += geometria.cbFromSegment(3, 0, 2)
cbe += geometria.cbFromSegment(1, 0, 1)
cbe += geometria.cbFromSegment(1, 0, 2)
geometria.cbe = cbe
O2 = PlaneStressNonLocal(geometria, E, v, b, l, z1,
                         Lr, af, fy=lambda x: -gamma*b)
O2.geometry.mask = None
O2.solveFromFile('U_2_Local.csv', plot=False)


_X, U1, U2, U3, UU = O.profile([0, h*0.49], [L, h*0.5], n=1000)
# plt.savefig('local.svg', transparent=True)
_X2, U12, U22, U32, UU2 = O2.profile([0, h*0.49], [L, h*0.5], n=1000)
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
# plt.plot(_X, w2(np.array(_X)), '-.', c='k', label='No Local Analítica')
plt.plot(_X2, np.array(UU2)[:, 1, 0], '--', c='k', label='Local')
plt.grid()
plt.xlabel('x')
plt.ylabel(r'$\Delta$')
plt.legend()
plt.show()

# plt.savefig('nolocal.svg', transparent=True)
