from FEM.EulerBernoulliBeam import EulerBernoulliBeam
from FEM.Geometry.Lineal import Lineal
import numpy as np
import matplotlib.pyplot as plt

L = 100
geometria = Lineal(L, 80, 1, 2)
geometria.cbe = [[0, 0.0], [1, 0.0], [-1, 0.0], [-2, 0.0]]

W = 10.0
Fx = 0.0
I = 1/12
A = 1
E = 30000000.0
EI = E*I
EA = E*A

O = EulerBernoulliBeam(geometria, EI, W)
O.solve()

# analitica = W*L**4/EI/8
print(O.U[-2][0])
plt.show()
a = 0
