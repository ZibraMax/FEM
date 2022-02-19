import numpy as np
import matplotlib.pyplot as plt
from FEM.EDO1D import EDO1D
from FEM.Geometry import Lineal


def a(x): return (x[0]**2-2)
def c(x): return (x[0]-3)
def f(x): return (x[0]**2-2)*6+(x[0]-3)*(3*x[0]**2)


cbe = [[0, 0], [-1, 3*1**2]]

lenght = 1
n = 500
o = 2
geometria = Lineal(lenght, n, o)
O = EDO1D(geometria, a, c, f)
O.cbe = cbe
O.solve()
plt.show()
