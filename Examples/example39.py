from numpy.testing._private.nosetester import get_package_name
from FEM.Torsion2D import Torsion2D
from FEM.Mesh.Delaunay import Delaunay
from FEM.Utils.polygonal import roundCorner, giveCoordsCircle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpltPath
import copy
r = 1.0
hueco = 0.4
N = 25
phi = 1.0
G = 80000.0

vert, seg1 = giveCoordsCircle([0, 0], r, n=N, sa=0, a=np.pi/2)
vert += [[0, r], [0, 0]]
a = seg1[-2][-1]
seg = copy.deepcopy(seg1)
seg[-1] = [a, a+1]
seg += [[a+1, 0]]
vertextra = [[hueco, 0], [hueco, hueco], [0, hueco]]
segextra = [[0, 1], [1, 2]]

holes = []
hole = {'center': [-2*r, 2*r], 'segments': segextra, 'vertices': vertextra}
holes += [hole]


params = Delaunay._strdelaunay(
    constrained=True, delaunay=True, a='0.005', o=2)

geometria = Delaunay(vert, params, nvn=1, holes_dict=holes)
GG = []
su = 0
for centroide in geometria.centroids:
    path = mpltPath.Path([[0, 0], [hueco, 0], [hueco, hueco], [0, hueco]])
    inside2 = path.contains_points([centroide])
    GG += [G]
    if inside2[0]:
        GG[-1] = 1
seg1[-1][-1] = a+1
geometria.segments = seg1
geometria.mask = None
O = Torsion2D(geometria, GG, phi)
O.geometry.show()
plt.show()
O.geometry.holes = None
O.solve()
plt.show()
integral = 0
for i, e in enumerate(O.elements):
    _, _u = e.giveSolution(domain='gauss-points')
    jac, dpz = e.J(e.Z.T)
    detjac = np.linalg.det(jac)
    integral += np.sum(_u*e.W*detjac)
print(integral*2/G)
