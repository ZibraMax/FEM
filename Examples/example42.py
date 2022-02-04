import numpy as np
import matplotlib.pyplot as plt
from FEM.NonLinealExample import NonLinealSimpleEquation
from FEM.Mesh.Lineal import Lineal

f0 = -1
Q = 0
u1 = 2**0.5


def a(x): return 1
def f(x): return f0


def u(x): return (1+x**2)**0.5
def dudx(x): return x/((1+x**2)**0.5)


cbn = [[0, 0]]
cbe = [[-1, u(1)]]

lenght = 1
n = 2
o = 1
geometria = Lineal(lenght, n, o)
O = NonLinealSimpleEquation(geometria, a, f)
O.cbe = cbe
O.cbn = cbn
O.solve()
x = np.linspace(0, 1, 100)
plt.gcf().get_axes()[0].plot(x, (1+x**2)**0.5)
plt.gcf().get_axes()[1].plot(x, 1/2*(1+x**2)**(-0.5)*2*x)
plt.show()
