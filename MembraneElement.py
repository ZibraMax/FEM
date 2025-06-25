import FEM
import numpy as np
import matplotlib.pyplot as plt

coords = np.array([[0.5, 0.0, 0.0],
                   [1.0, 1.0, 1.0],
                   [0.5, 1.0, 1.0],
                   [0.0, 0.5, 0.0]])

gdl = np.array([[0, 1, 2],
                [3, 4, 5],
                [6, 7, 8],
                [9, 10, 11]]).T
m = gdl.shape[0]*gdl.shape[1]
E = 210e3
nu = 0.3

c11 = E / (1 - nu**2)  # Revisar
c12 = nu * c11  # Revisar
c22 = c11
c66 = E / (2 * (1 + nu))
C = np.array([
    [c11, c12, 0.0],
    [c12, c22, 0.0],
    [0.0, 0.0, c66]])
t = 0.05
e = FEM.Quadrilateral(coords, gdl, 3, boundary=True)
__x, __p = e.T(e.domain.T)
_x, _p = e.T(e.Z.T)
_j, _dp = e.J(e.Z.T)
weights = e.W
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

Ke = np.zeros((m, m))

ax.plot_trisurf(*__x.T, alpha=0.3)
for x, jac, wi, ni, dni in zip(_x, _j, weights, _p, _dp):
    JP = np.linalg.pinv(jac)
    e3 = np.cross(*jac, axis=0)
    detjac = np.linalg.norm(e3)
    e1 = jac[0]
    e2 = jac[1]
    e1 /= np.linalg.norm(e1)
    e2 /= np.linalg.norm(e2)
    e3 /= np.linalg.norm(e3)
    dpx = JP @ dni  # Oloverga!
    _m = len(e.gdl.T)
    o = [0.0]*_m
    B = np.array([
        [*dpx[0, :], *o],
        [*o, *dpx[1, :]],
        [*dpx[1, :], *dpx[0, :]]])

    Ke += t*(B.T@C@B)*detjac*wi

    ax.plot(*x, 'o', c='k', zorder=200)
    ax.quiver(*x, *e3/3, color="green", zorder=200)
    ax.quiver(*x, *e1/3, color="blue", zorder=200)
    ax.quiver(*x, *e2/3, color="red", zorder=200)
plt.show()
