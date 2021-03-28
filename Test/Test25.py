from FEM.Mesh.Geometry import Geometry
from FEM.Mesh.Delaunay import Delaunay
from FEM.PlaneStrain import PlaneStrain
from FEM.Utils.polygonal import roundCorner, giveCoordsCircle
import matplotlib.pyplot as plt
import numpy as np

E = 30*10**(5)
v = 0.25

b = 10
h = 20
he = h/4
ancho_en_h10_in = 18
ancho_en_h20_in = 10
p0 = 200
pp = 1000
ppx = pp*3/5
ppy = -pp*4/5


def darPolinomio(X, Y):
    n = len(X)
    A = np.zeros([n, n])
    B = np.zeros([n, 1])
    for i in range(n):
        for j in range(n):
            A[i, j] = X[i]**j
        B[i, 0] = Y[i]

    U = np.linalg.solve(A, B)

    def f(x):
        suma = 0
        for i in range(n):
            suma += U[i, 0]*x**i
        return suma
    return f


n = 20

parabola = darPolinomio(np.array([0, 10, 20]), np.array(
    [0, b-ancho_en_h10_in/2, b-ancho_en_h20_in/2]))

c = [
    [0, 0],
    [2*b, 0]]

for i in range(1, n):
    x = 2*b-parabola(h/n*i)
    y = h/n*i
    c += [[x, y]]

c += [[2*b-parabola(4*he), 4*he],
      [parabola(4*he), 4*he]]
for i in reversed(range(1, n)):
    x = parabola(h/n*i)
    y = h/n*i
    c += [[x, y]]

holes = []
radi = 2
cent = [b, h/2]
vert, seg = giveCoordsCircle(cent, radi, n=50)
hole = {'center': cent, 'segments': seg, 'vertices': vert}
holes += [hole]
params = Delaunay._strdelaunay(constrained=True, delaunay=True, a='0.1', o=2)
geometria = Delaunay(c, params, nvn=2, holes_dict=holes)
geometria.generateSegmentsFromCoords([0, 0], [2*b, 0])
geometria.generateSegmentsFromCoords(
    [2*b-parabola(4*he), 4*he], [parabola(4*he), 4*he])
geometria.cbe = geometria.cbFromSegment(-2, 0, 1)
geometria.cbe += geometria.cbFromSegment(-2, 0, 2)
geometria.saveMesh('Mesh_tests/tunel')
geometria.show()
plt.show()
geometria.loadOnSegment(-1, fy=lambda s: -p0)
geometria.mask = None
O = PlaneStrain(geometria, E, v)
O.elementMatrices()
O.ensembling()
O.borderConditions()
O.solveES()
O.postProcess()

plt.show()
