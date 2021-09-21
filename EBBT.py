from FEM.EulerBernoulliBeam import EulerBernoulliBeam
from FEM.Mesh.Lineal import Lineal
import numpy as np
import matplotlib.pyplot as plt

L = 3
geometria = Lineal(L, 100, 1, 2)
geometria.cbe = [[0, 0.0], [1, 0.0]]

W = 10.0
EI = 100.0

O = EulerBernoulliBeam(geometria, EI, W)
O.solve()

analitica = W*L**4/EI/8
print(analitica, O.U[-2])

plt.show()
