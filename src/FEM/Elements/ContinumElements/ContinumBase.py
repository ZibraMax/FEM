from ..Element import Element
import numpy as np


class ContinumBase(Element):
    """Base class for continuum elements.

    Args:
        coords (np.ndarray): Element coordinate matrix
        _coords (np.ndarray): Element coordinate matrix for graphical interface purposes
        gdl (np.ndarray): Degree of freedom matrix
    """

    def __init__(self, coords: np.ndarray, _coords: np.ndarray, gdl: np.ndarray, **kargs) -> None:
        """Create a continuum element

        Args:
            coords (np.ndarray): Element coordinate matrix
            _coords (np.ndarray): Element coordinate matrix for graphical interface purposes
            gdl (np.ndarray): Degree of freedom matrix
        """
        super().__init__(coords, _coords, gdl, **kargs)

        # Initialize element matrices and vectors
        self.Ke = None  # Stiffness matrix
        self.Fe = None  # Force vector
        self.Ue = None  # Displacement vector
        self.Qe = None  # Internal force vector

    def elementMatrices(self) -> None:
        """Calculate element matrices and vectors.

        This method should be implemented in derived classes.
        """
        raise NotImplementedError(
            "This method should be implemented in derived classes.")
