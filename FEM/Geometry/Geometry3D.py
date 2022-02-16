from .Geometry import Geometry
from typing import Callable
import numpy as np

class Geometry3D(Geometry):
    """Creates a 3D geometry"""
    def __init__(self, dictionary: list, gdls: list, types: list, nvn: int = 1, segments: list = [], fast=False):
        Geometry.__init__(self, dictionary, gdls, types, nvn, segments, fast)

    def show(self, texto: int = 10, bolita: int = 0, draw_segs: bool = True, draw_labels: bool = False, draw_bc: bool = False, label_bc: bool = False) -> None:
        pass