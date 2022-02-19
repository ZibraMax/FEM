from FEM.Elements import E1D
from FEM.Elasticity2D import PlaneStressOrthotropic as pso
from FEM.Geometry.Geometry import Geometry
import matplotlib.pyplot as plt

a = 120.0
b = 160.0
p0 = 10.0

E1 = 31*10**6
E2 = 2.7*10**6
G12 = 0.75*10**6
v12 = 0.28
t = 0.036

gdls = [[0.0, 0.0], [a, 0.0], [a, b], [0, b]]
diccs = [[0, 1, 2], [0, 2, 3]]
types = ['T1V', 'T1V']
nvn = 2
regions = [[0, 1], [1, 2], [2, 3], [3, 0]]
geometry = Geometry(dictionary=diccs, gdls=gdls,
                    types=types, nvn=nvn, regions=regions)

cb = geometry.cbFromRegion(3, 0.0, 1)
cb += geometry.cbFromRegion(3, 0.0, 2)
geometry.setCbe(cb)
def f(s): return p0


geometry.loadOnRegion(1, f)
# geometry.show(draw_labels=True, draw_bc=True, label_bc=True)
# plt.show()
OFEM = pso(geometry, E1, E2, G12, v12, t)
OFEM.solve()
plt.show()
