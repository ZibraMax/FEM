from FEM.Elements.E3D.Brick import Brick
import numpy as np
coords = np.array([[0, 0, 0.0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [
                  0, 0, 1.0], [1, 0, 1], [1.5, 1.5, 1.5], [0, 1, 1]])
coords2 = coords + 1.5
gdl = np.array([[0, 0, 0, 0, 0, 0, 0, 0]])
e = Brick(coords=coords, gdl=gdl)
e2 = Brick(coords=coords2, gdl=gdl)
domain = e.T(e.domain.T)[0]
domain2 = e2.T(e.domain.T)[0]
points = np.array(
    np.array([[-1, 0, 0], [0.5, 0.5, 0.5]]).tolist()+domain.tolist())

r1 = e.isInside(points)
r2 = e.isInside(domain2)
a = 0
