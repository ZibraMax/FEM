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
            E = np.array([E]*len(geometry.elements))
        if isinstance(A, float) or isinstance(A, int):
            A = np.array([A]*len(geometry.elements))

        self.properties["E"] = E
        self.properties["A"] = A
        self.ngdl = len(self.geometry.elements)
        self.A = sparse.lil_matrix((self.ngdl, self.ngdl))
        self.G = np.zeros(self.ngdl)

    def elementMatrices(self):
        """Calculate element stiffness matrices and assemble global stiffness matrix.
        """
        for i, e in enumerate(self.geometry.elements):
            # Element length and angle
            x1 = self.geometry.nodes[e[0]]  # 2D or 3D coordinates
            x2 = self.geometry.nodes[e[1]]  # 2D or 3D coordinates
            vector = x2 - x1
            basis = np.array([1, 0, 0])
            angle = np.atan2(np.linalg.norm(
                np.cross(vector, basis)), np.dot(vector, basis))
