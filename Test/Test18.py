import matplotlib.pyplot as plt
import numpy as np

from FEM.Elements.E2D.LTriangular import LTriangular
from FEM.Elements.E2D.QTriangular import QTriangular
from FEM.Elements.E2D.Serendipity import Serendipity
from FEM.Elements.E2D.Quadrilateral import Quadrilateral

coords = np.array([[0.5, 0.5], [0.8, 0.5], [0.65, 0.8]])
gdls = np.array([[0, 1, 2]])
e = LTriangular(coords, gdls)


coords = coords.tolist()
for i in range(len(coords)-1):
    coords += [[(coords[i][0]+coords[i+1][0])/2,
                (coords[i][1]+coords[i+1][1])/2]]
coords += [[(coords[i+1][0]+coords[0][0])/2, (coords[i+1][1]+coords[0][1])/2]]
gdl = np.array([[1, 2, 3, 4, 5, 6]])
e2 = QTriangular(coords, gdl)

coords = [[1, 1], [3, 2], [3.5, 3], [0.5, 4], [2, 1.5],
          [3.25, 2.5], [(3.5+.5)/2, 3.5], [(0.5+1)/2, 5/2]]
gdl = np.array([[1, 2, 3, 4, 5, 6, 7, 8]])
e3 = Serendipity(coords, gdl)

coords = [[1, 1], [3, 2], [3.5, 3], [0.5, 4]]
gdl = np.array([[1, 2, 3, 4]])
e4 = Quadrilateral(coords, gdl)

e3.borders[0].draw()
plt.show()
