import numpy as np
import matplotlib.pyplot as plt
from FEM.Geometry.Geometry import Geometry
from FEM.Poisson2D import Poisson2D
phi = 0
nodos = [[0.5, 0.0], [1.0, 0.0], [0.5, 0.5],
         [1.0, 0.5], [0.5, 1.0], [1.0, 1.0]]
elementos = [[0, 1, 3, 2], [2, 3, 5, 4]]
tipos = ['C1V', 'C1V']
geometry = Geometry(elementos, nodos, tipos)
geometry.cbe = [[0, 1], [1, 0], [3, 1/4], [5, 1], [4, 1]]
geometry.cbn = [[2, 0]]
# geometry.show(draw_segs=True, draw_labels=True, draw_bc=True, label_bc=True)
# plt.show()
O = Poisson2D(geometry, phi)
O.solve(plot=False)
# plt.show()
print(O.U[2, 0])
plt.show()


elementos = [[0, 1, 3], [0, 3, 2], [2, 3, 5], [2, 5, 4]]
tipos = ['T1V', 'T1V', 'T1V', 'T1V']
geometry = Geometry(elementos, nodos, tipos)
geometry.cbe = [[0, 1], [1, 0], [3, 1/4], [5, 1], [4, 1]]
geometry.cbn = [[2, 0]]
geometry.show(draw_segs=True, draw_labels=True, draw_bc=True, label_bc=True)
plt.show()
O = Poisson2D(geometry, phi)
O.solve()
print(O.U[2, 0])
plt.show()
elementos = [[0, 1, 2], [2, 1, 3], [2, 3, 5], [2, 5, 4]]
tipos = ['T1V', 'T1V', 'T1V', 'T1V']
geometry = Geometry(elementos, nodos, tipos)
geometry.cbe = [[0, 1], [1, 0], [3, 1/4], [5, 1], [4, 1]]
geometry.cbn = [[2, 0]]
geometry.show(draw_segs=True, draw_labels=True, draw_bc=True, label_bc=True)
plt.show()
O = Poisson2D(geometry, phi)
O.solve()
print(O.U[2, 0])
plt.show()
