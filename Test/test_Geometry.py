"""Geometry tests"""

from FEM.Utils import polygonal
from FEM.Geometry.Geometry import Geometry
import unittest
import os
import sys
import inspect
FILENAME = 'Test/resources/beam.msh'


class TestGeometry(unittest.TestCase):
    """Test the geometry class"""

    def test_meshRect(self):
        """Creation a rectangular geometry file
        """
        h = 0.5
        L = 2.0
        nex = 40
        ney = 10
        polygonal.enmalladoFernando(L, h, nex, ney, FILENAME)
        geometry = Geometry.loadmsh(FILENAME)
        self.assertEqual(len(geometry.elements), nex*ney)
        self.assertEqual(len(geometry.segments), 0)

    def test_generateSegmentFromCoords(self):
        """Test if the generated geometry can be modified with segments from coordinates. Then saves the new geometry to a file
        """
        h = 0.5
        L = 2.0
        geometry = Geometry.loadmsh(FILENAME)
        geometry.generateSegmentsFromCoords([0, 0], [L, 0])
        geometry.generateSegmentsFromCoords([L, 0], [L, h])
        geometry.generateSegmentsFromCoords([L, h], [0, h])
        geometry.generateSegmentsFromCoords([0, h], [0, 0])
        self.assertEqual(len(geometry.segments), 4)
        geometry.saveMesh('Test/resources/beam_ws')


if __name__ == '__main__':
    currentdir = os.path.dirname(os.path.abspath(
        inspect.getfile(inspect.currentframe())))
    parentdir = os.path.dirname(currentdir)
    sys.path.insert(0, parentdir)
    from FEM.Geometry.Geometry import Geometry
    unittest.main()
