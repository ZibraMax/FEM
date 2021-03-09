# Creación de elementos
import matplotlib.pyplot as plt
import numpy as np

import FEM

ELEMENTOS = []
coords = [[1,1],[3,2],[3.5,3],[0.5,4],[2,1.5],[3.25,2.5],[(3.5+.5)/2,3.5],[(0.5+1)/2,5/2]]
gdl = np.array([[1,2,3,4,5,6,7,8]])
eRS = FEM.Elements.E2D.Serendipity(coords,gdl)
ELEMENTOS.append(eRS)

coords = [[1,1],[3,2],[3.5,3],[0.5,4]]
gdl = np.array([[1,2,3,4]])
eRC = FEM.Elements.E2D.Quadrilateral(coords,gdl)
ELEMENTOS.append(eRC)
coords = [[-1.5,-1],[1,0],[0,5]] #[[-1.5, -1], [1, 0], [0, 5], [-0.25, -0.5], [0.5, 2.5], [-0.75, 2.0]]
for i in range(len(coords)-1):
	coords+=[[(coords[i][0]+coords[i+1][0])/2,(coords[i][1]+coords[i+1][1])/2]]
coords+=[[(coords[i+1][0]+coords[0][0])/2,(coords[i+1][1]+coords[0][1])/2]]
gdl = np.array([[1,2,3,4,5,6]])
eTC = FEM.Elements.E2D.QTriangular(coords,gdl)
ELEMENTOS.append(eTC)

coords = [[1.4,1.5],[3.1,2.4],[2,3]]
gdl = np.array([[1,2,3]])
eTL = FEM.Elements.E2D.LTriangular(coords,gdl)

coords = [[3],[4],[5]]
gdl = np.array([[3,4,5]])
e = FEM.Elements.E1D.QuadraticElement(coords,gdl,3)
ELEMENTOS.append(e)
coords = [[3],[5]]
gdl = np.array([[3,5]])
e = FEM.Elements.E1D.LinealElement(coords,gdl,3)
ELEMENTOS.append(e)

# Dibujar elementos
# for e in ELEMENTOS:
# 	e.draw()
# 	plt.show()


# Esta Adentro
# print(e.isInside(np.array([1,2,3,4,5])))


# Inverse Mapping
z = eTL.inverseMapping(np.array([[1.3,2.5,3.5],[1.5,2.6,8.5]]))
# z = eTL.inverseMapping(np.array([[1.3],[1.5]]))
print(z)
# print(eTL.isInside(np.array([3.5,2.5])))


# Jacobian Graph
# for e in ELEMENTOS:
# 	e.jacobianGraph()
# 	plt.show()

