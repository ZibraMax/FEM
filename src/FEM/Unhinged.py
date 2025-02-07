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
            x1, y1 = self.geometry.nodes[e[0]]
            x2, y2 = self.geometry.nodes[e[1]]
            L = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            c = (x2 - x1) / L
            s = (y2 - y1) / L

            # Element stiffness matrix
            k = self.properties["E"][i] * self.properties["A"][i] / L
            k_matrix = k * np.array([[c**2, c*s, -c**2, -c*s],
                                     [c*s, s**2, -c*s, -s**2],
                                     [-c**2, -c*s, c**2, c*s],
                                     [-c*s, -s**2, c*s, s**2]])

            # Assemble global stiffness matrix
            self.A[i:i+4, i:i+4] += k_matrix
