import numpy as np
import matplotlib.pyplot as plt
import FEM

ELEMENTOS = []
coords = [[1,1],[3,2],[3.5,3],[0.5,4],[2,1.5],[3.25,2.5],[(3.5+.5)/2,3.5],[(0.5+1)/2,5/2]]
gdl = [1,2,3,4,5,6,7,8]
eRS = FEM.Elements.E2D.Serendipity(coords,gdl)
ELEMENTOS.append(eRS)

coords = [[1,1],[3,2],[3.5,3],[0.5,4]]
gdl = [1,2,3,4]
eRC = FEM.Elements.E2D.Quadrilateral(coords,gdl)
ELEMENTOS.append(eRC)
coords = [[-1.5,-1],[1,0],[0,5]]
for i in range(len(coords)-1):
	coords+=[[(coords[i][0]+coords[i+1][0])/2,(coords[i][1]+coords[i+1][1])/2]]
coords+=[[(coords[i+1][0]+coords[0][0])/2,(coords[i+1][1]+coords[0][1])/2]]
gdl = [1,2,3,4,5,6]
eTC = FEM.Elements.E2D.QTriangular(coords,gdl)
ELEMENTOS.append(eTC)

coords = [[1.4,1.5],[3.1,2.4],[2,3]]
gdl = [1,2,3]
eTL = FEM.Elements.E2D.LTriangular(coords,gdl)
ELEMENTOS.append(eTL)


for e in ELEMENTOS:
	e.draw()
	plt.show()
	e.jacobianGraph()
	plt.show()


# coords = [[3],[4],[5]]
# gdl = [3,4,5]
# e = FEM.Elements.E1D.QuadraticElement(coords,gdl)
# print(*e.T(np.array([-1,0,1])))
# print(e.dpsis(np.array([[-1,0],[0,0]])))
# print(e.dpsis([-1,0]))
# print(e.dpsis([0,0]))

# print(eRS.J(np.array([[0,1],[0,0]])))
# print(e.J(np.array([i*0.1 for i in range(10)])))
# e.draw()
# plt.show()
