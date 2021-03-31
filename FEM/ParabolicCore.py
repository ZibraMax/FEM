"""Transitory analysis core
"""
from .Core import Core, np


class ParabolicCore(Core):
    def __init__(self, geometry):
        Core.__init__(self, geometry)
        self.M = np.zeros([self.ngdl, self.ngdl])
        self.F = np.zeros([self.ngdl, 1])
