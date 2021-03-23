"""Transitory analysis core
"""
from .Core import Core


class Transitory(Core):
    def __init__(self, geometry):
        Core.__init__(self, geometry)
