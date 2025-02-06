"""Origami Bar and Hinge Model
"""

from .Core import Core, Geometry, logging
from .Solvers import LinealSparse
from typing import Callable, Tuple
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from scipy import sparse  # Creation of sparse matrices


class Truss(Core):
    """Base truss analysis using the Linear force method. 
    """

    def __init__(self, geometry: Geometry, E: float, A: float, solver: LinealSparse = LinealSparse(), **kwargs):
        """Initialize the truss model.

        Args:
            geometry (Geometry): Geometry object containing nodes and elements.
            solver (LinealSparse): Solver object for linear analysis.
        """
        super().__init__(geometry, solver, **kwargs)
        self.name = 'Truss analysis (2D and 3D)'
        if isinstance(E, float) or isinstance(E, int):
            E = [E]*len(geometry.elements)
        if isinstance(A, float) or isinstance(A, int):
            A = [A]*len(geometry.elements)

        self.properties["E"] = E
        self.properties["A"] = A
