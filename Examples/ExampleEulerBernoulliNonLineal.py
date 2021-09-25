from FEM.EulerBernoulliBeam import EulerBernoulliBeamNonLineal
from FEM.Mesh.Lineal import Lineal
import numpy as np
import matplotlib.pyplot as plt

L = 100
geometria = Lineal(L, 80, 1, 3)
geometria.cbe = [[0, 0.0], [1, 0.0], [2, 0.0], [-1, 0.0], [-2, 0.0], [-3, 0.0]]

W = 10.0
Fx = 0.0
I = 1/12
A = 1
E = 30000000.0
EI = E*I
EA = E*A

O = EulerBernoulliBeamNonLineal(geometria, EI, EA, Fx, W)
O.solve()

# analitica = W*L**4/EI/8
print(O.U[-2][0])
plt.show()
a = 0
