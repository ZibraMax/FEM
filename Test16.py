import numpy as np
import matplotlib.pyplot as plt
import FEM
from FEM import Mesh

E = 2.1*10**6
v = 0.2
u = 0.001
a = 5
t = 0.5
l = 0.1
z1 = 0.5
geometria = Mesh.Geometry.loadmsh('Mesh_tests/EnmalladoTesis.msh')
geometria.generateSegmentsFromCoords([0,0],[a,0])
geometria.generateSegmentsFromCoords([0,a],[a,a])
geometria.generateSegmentsFromCoords([a,a],[0,a])
geometria.generateSegmentsFromCoords([0,a],[0,0])
geometria.cbFromSegment(3,0,1)
geometria.cbFromSegment(3,0,2)
geometria.cbFromSegment(1,u,1)
def af(l0,rho):
	return l0*np.exp(-rho)
O = FEM.PlaneStressNonLocal(geometria,E,v,t,l,z1,Lr=6*l,af=af)
O.solve()
plt.show()