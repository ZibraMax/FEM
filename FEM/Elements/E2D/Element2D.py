from ..Element import *
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.path as mpltPath


class Element2D(Element):
    def __init__(self, coords, _coords, gdl):
        Element.__init__(self, coords, _coords, gdl)
        self._coordsg = np.array(
            self._coords.tolist()+[self._coords[0].tolist()])

    def draw(self):
        _z = self.domain
        _x, _p = self.T(_z.T)
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        l = []
        l.append('Element')
        l.append('Nodes')
        for i in range(self.n):
            surf = ax.plot_trisurf(*_x.T, _p[:, i], alpha=0.3)
            surf._facecolors2d = surf._facecolors3d
            surf._edgecolors2d = surf._edgecolors3d
            l.append(r'$\psi_{'+format(i)+r'}$')
        __coords = np.array(self._coords.tolist()+[self._coords[0].tolist()]).T
        ax.plot(*__coords, [0]*len(__coords.T), '-', color='black')
        ax.plot(*self.coords.T, [0]*len(self.coords), 'o', color='blue')
        ax.legend(l)

    def jacobianGraph(self):
        _z = self.domain
        _x, _p = self.T(_z.T)
        _j = self.J(_z.T)[0]
        __j = np.linalg.det(_j)
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        l = []
        surf = ax.plot_trisurf(*_x.T, __j, cmap='magma')
        surf._facecolors2d = surf._facecolors3d
        surf._edgecolors2d = surf._edgecolors3d
        l.append('Element')
        l.append('Nodes')
        l.append(r'$|J|$')
        cbar = fig.colorbar(surf)
        __coords = np.array(self._coords.tolist()+[self._coords[0].tolist()]).T
        ax.plot(*__coords, [0]*len(__coords.T), '-', color='black')
        ax.plot(*self.coords.T, [0]*len(self.coords), 'o', color='blue')
        ax.legend(l)

    def isInside(self, x):  # TODO hacer que esto sea vectorizado. En teoría ya lo es xd
        path = mpltPath.Path(self._coords)
        inside2 = path.contains_points([x])
        return inside2[0]
