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

    def set_constitutive_model(self, CM) -> None:
        self.constitutive_model = CM

    def calculate_deformation_gradient(self, dpt):
        t_0 = self.coords.copy()
        t_t = t_0 + self.Ue.T
        F = dpt @ t_t
        return F

    def green_lagrange_strain(self, F):
        """Calculate the Green-Lagrange strain tensor from the deformation gradient."""
        E = 0.5 * (F.T @ F - np.eye(F.shape[0]))
        return E

    def calculate_BNL(self, dpx) -> np.ndarray:
        """Calculate the BNL matrix for the element.

        This method should be implemented in derived classes.
        """
        raise NotImplementedError(
            "This method should be implemented in derived classes.")

    def calculate_BL(self, dpx) -> np.ndarray:
        """Calculate the BL matrix for the element.

        This method should be implemented in derived classes.
        """
        raise NotImplementedError(
            "This method should be implemented in derived classes.")

    def organize_S(self, S) -> tuple:
        """Organize the stress tensor S and the constitutive matrix C.

        This method should be implemented in derived classes.
        """
        raise NotImplementedError(
            "This method should be implemented in derived classes.")

    def transformation_matrix(self) -> np.ndarray:
        """Calculate the transformation matrix for the element.

        This method should be implemented in derived classes.
        """
        raise NotImplementedError(
            "This method should be implemented in derived classes.")

    def get_local_jacobian(self, jac: np.ndarray) -> np.ndarray:
        """Get the local Jacobian matrix for the element.

        Args:
            jac (np.ndarray): The Jacobian matrix from the integration points.

        Returns:
            np.ndarray: The local Jacobian matrix.
        """
        raise NotImplementedError(
            "This method should be implemented in derived classes.")

    def elementMatrices(self) -> None:
        """Calculate element matrices and vectors.

        This method should be implemented in derived classes.
        """

        """Calculate the lineal stiffness matrix.

        This method should be implemented in derived classes.
        """
        _x = self._x
        _p = self._p
        _j = self._j
        _dp = self._dp
        weights = self.W

        K = 0.0
        F = 0.0
        for x, jac, wi, ni, dni in zip(_x, _j, weights, _p, _dp):
            J = self.get_local_jacobian(jac)
            detjac = np.linalg.det(J)
            # Esto asume que el jacobiano tiene inversa
            dpx = np.linalg.inv(J) @ dni
            # du = self.Ue.T @ dpx  # Creo
            F = self.calculate_deformation_gradient(dpx)
            E = self.green_lagrange_strain(F)
            C, S = self.constitutive_model(E)

            S_stiff, S_force = self.organize_S(S)

            BL = self.calculate_BL(dpx)
            BNL = self.calculate_BNL(dpx)
            K += (BL.T @ C @ BL + BNL.T @ S_stiff @ BNL) * detjac * wi
            # De donde sale ese du? ->
            F += (BL.T @ S_force) * detjac * wi

        T = self.transformation_matrix()
        if T is not None:
            K = T.T @ K @ T
            F = T.T @ F
        return K, F
