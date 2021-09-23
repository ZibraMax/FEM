"""
***************************************
Numerical validation of the ODE2D Class
***************************************

Tests
#####

"""

from FEM.Mesh.Geometry import Geometry
from FEM.ODE2D import ODE2D
import matplotlib.pyplot as plt
import unittest
import os
import numpy as np
import sys
import inspect
TOL = 0.1


class TestODE2D(unittest.TestCase):
    """Parial differential equations in 2D tests
    """

    def test_polinomyal_function(self):
        pass


if __name__ == '__main__':
    currentdir = os.path.dirname(os.path.abspath(
        inspect.getfile(inspect.currentframe())))
    parentdir = os.path.dirname(currentdir)
    sys.path.insert(0, parentdir)
    from FEM.Mesh.Geometry import Geometry
    from FEM.EDO2D import EDO2D
    unittest.main()
