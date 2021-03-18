from ..Element import *
from .LinearScheme import *
import numpy as np
import matplotlib.pyplot as plt


class Element1D(Element, LinearScheme):
    def __init__(self, coords, gdl, n):
        coords = np.array(coords).reshape([len(coords), 1])
        _coords = np.array([coords[0][0], coords[-1][0]])
        Element.__init__(self, coords, _coords, gdl)
        LinearScheme.__init__(self, n)

    def draw(self, m=100):
        _z = self.domain
        _x, _p = self.T(_z)
        _y = 0
        for i in range(self.n):
            plt.plot(_x, _p[:, i], '--', label=r'$\psi_{'+format(i)+r'}$')
            _y += _p[:, i]
        plt.plot(_x, [0]*len(_x), '-', color='black', label='Element')
        plt.plot(self.coords.T[0], [0]*len(self.coords.T[0]),
                 'o', color='blue', label='Nodes')
        plt.plot(_x, _y, '-.', label=r'$\sum_{i=0}^{n}\psi$')
        plt.legend()
        plt.grid()

    def jacobianGraph(self):
        pass  # No implementado porque el jacobiano es constante en cualquier elemento lienal

    def isInside(self, x):
        return (x >= self.coords.T[0][0])*(x <= self.coords.T[0][-1])
