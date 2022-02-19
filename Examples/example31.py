from FEM.Geometry.Geometry import Geometry
from FEM.Heat2D import Heat2D
import matplotlib.pyplot as plt
b = 4
h = 2

kx = 10
ky = 15
Ta = 100
beta = 5
T0 = 0

nx = 4
ny = 2
hx = b/nx
hy = h/ny
gdls = []

for j in range(ny+1):
    for i in range(nx+1):
        gdls += [[i*hx, j*hy]]
elementos = []
for j in range(ny):
    for i in range(nx):  # i = 0 j = 0
        elementos += [[i+j*(nx+1), i+1+j*(nx+1), i+1+(j+1)*(nx+1), i+(j+1)
                       * (nx+1)]]
tipos = ['C1V']*len(elementos)
geometria = Geometry(elementos, gdls, tipos, 1)
geometria.generateRegionFromCoords([0, 0], [b, 0])
geometria.generateRegionFromCoords([b, 0], [b, h])
geometria.generateRegionFromCoords([b, h], [0, h])
geometria.generateRegionFromCoords([0, h], [0, 0])
cbe = geometria.cbFromSegment(2, T0)
geometria.setCbe(cbe)

O = Heat2D(geometria, kx, ky)
O.defineConvectiveBoderConditions(1, beta, Ta)
O.geometry.show(draw_segs=True, draw_bc=True)
plt.show()
O.solve()
plt.show()
O.print()
