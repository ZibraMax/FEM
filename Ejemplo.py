from FEM.EDO1D import EDO1D
from FEM.Mesh.Lineal import Lineal
import numpy as np
import matplotlib.pyplot as plt

L = 10

P = 30000
E = 20000000
A0 = 4.5
A1 = 3*A0
def EA(x): return E*(A0+x/L*(A1-A0))


def sf(x): return 0


n = 5
geometry = Lineal(L, n, 1)
geometry.cbe = [[-1, 0]]
geometry.cbn = [[0, P]]
O = EDO1D(geometry, EA, lambda x: 0, sf)
O.solve()
plt.show()
