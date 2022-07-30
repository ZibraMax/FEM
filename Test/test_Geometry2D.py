"""Geometry tests"""
from FEM.Geometry.Geometry import Geometry2D
from FEM.Geometry import Region1D
import unittest
FILENAME = 'Test/resources/beam.json'


class TestGeometry2D(unittest.TestCase):
    """Test the geometry class"""

    def test_importJSON(self):
        """Import of geometry using file
        """
        geometry = Geometry2D.importJSON(FILENAME)
        self.assertEqual(len(geometry.regions), 0)

    def test_addRegions(self):
        """Test if the generated geometry can be modified with regions from coordinates. Then saves the new geometry to a file
        """
        h = 0.5
        L = 2.0
        regions = []
        regions.append(Region1D([[0, 0], [L, 0]]))
        regions.append(Region1D([[L, 0], [L, h]]))
        regions.append(Region1D([[L, h], [0, h]]))
        regions.append(Region1D([[0, h], [0, 0]]))
        geometry = Geometry2D.importJSON(FILENAME)
        geometry.addRegions(regions)
        self.assertEqual(len(geometry.regions), 4)
        geometry.exportJSON('Test/resources/beam_ws.json')

