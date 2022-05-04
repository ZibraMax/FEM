"""Geometry tests"""

from FEM.Utils import polygonal
from FEM.Geometry.Geometry import Geometry
import unittest
FILENAME = 'Test/resources/beam.json'


class TestGeometry(unittest.TestCase):
    """Test the geometry class"""

    def test_meshRect(self):
        """Creation a rectangular geometry file
        """
        h = 0.5
        L = 2.0
        nex = 40
        ney = 10
        coords, dicc = polygonal.enmalladoFernando(L, h, nex, ney)
        geometry = Geometry(dicc, coords, ['C2V']*len(dicc), nvn=2)
        geometry.exportJSON(FILENAME)
        self.assertEqual(len(geometry.elements), nex*ney)
        self.assertEqual(len(geometry.regions), 0)


if __name__ == '__main__':
    unittest.main()
