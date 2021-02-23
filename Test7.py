import numpy as np
import matplotlib.pyplot as plt
import FEM
from FEM import Mesh

a = lambda x: (x**2-2)
c = lambda x: (x-3)
f = lambda x: (x**2-2)*6+(x-3)*(3*x**2)
cbe = [[0,0],[-1,3*1.5**2]]

lenght = 1.5
n = 10
o = 2

geometria = Mesh.Lineal(lenght, n, o)
O = FEM.EDO1D(geometria,a,c,f)
O.cbe = cbe
O.elementMatrices()
O.ensembling()
O.borderConditions()
O.solveES()
O.postProcess()
plt.show()
