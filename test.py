import numpy as np
import matplotlib.pyplot as plt
import FEM

ELEMENTOS = []
coords = [[1,1],[3,2],[3.5,3],[0.5,4],[2,1.5],[3.25,2.5],[(3.5+.5)/2,3.5],[(0.5+1)/2,5/2]]
gdl = [[1,2,3,4,5,6,7,8]]
eRS = FEM.Elements.E2D.Serendipity(coords,gdl)
ELEMENTOS.append(eRS)

coords = [[1,1],[3,2],[3.5,3],[0.5,4]]
gdl = [[1,2,3,4]]
eRC = FEM.Elements.E2D.Quadrilateral(coords,gdl)
ELEMENTOS.append(eRC)
coords = [[-1.5,-1],[1,0],[0,5]]
for i in range(len(coords)-1):
	coords+=[[(coords[i][0]+coords[i+1][0])/2,(coords[i][1]+coords[i+1][1])/2]]
coords+=[[(coords[i+1][0]+coords[0][0])/2,(coords[i+1][1]+coords[0][1])/2]]
gdl = [[1,2,3,4,5,6]]
eTC = FEM.Elements.E2D.QTriangular(coords,gdl)
ELEMENTOS.append(eTC)

coords = [[1.4,1.5],[3.1,2.4],[2,3]]
gdl = [[1,2,3]]
eTL = FEM.Elements.E2D.LTriangular(coords,gdl)
ELEMENTOS.append(eTL)
Ue = np.array([[2,3,4],[1,2,3],[4,5,6]])
eTL.Ue = Ue
print(eTL.Z.T)
print(eTL.J(eTL.Z.T)[0])


coords = [[3],[4],[5]]
gdl = [[3,4,5]]
e = FEM.Elements.E1D.QuadraticElement(coords,gdl,3)
# print(e.isInside(np.array([1,2,3,4,5])))
# Ue = np.array([[2,3,4],[1,2,3],[4,5,6]])
# e.Ue = Ue
# print(e.giveSolution())
# print(eRC.T(np.array([[1.3,1.2],[2.67,3]]))[0])
# z = eRC.inverseMapping(np.array([[1.3,1.2],[2.67,3]]))
# z = eRC.inverseMapping([1.3,1.2])
# print(eRC.T(z)[0])

# for e in ELEMENTOS:
# 	e.draw()
# 	plt.show()
# 	e.jacobianGraph()
# 	plt.show()


# print(*e.T(np.array([-1,0,1])))
# print(e.dpsis(np.array([[-1,0],[0,0]])))
# print(e.dpsis([-1,0]))
# print(e.dpsis([0,0]))

# print(eRS.J(np.array([[0,1],[0,0]])))
# print(e.J(np.array([i*0.1 for i in range(10)])))
# z = e.inverseMapping(3.5)
# print(z)
# print(e.T(z)[0])

#TODO revisar inverse mapping múultiple elementos lineales porque hay que meterlos como vectores verticales
# z = e.inverseMapping(np.array([[3.5],[4.5]]))
# print(z)
# print(e.T(z)[0])
# e.draw()
# plt.show()
