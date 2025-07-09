import numpy as np


class ContinumBase():
    """Base class for continuum elements.

    Args:
        coords (np.ndarray): Element coordinate matrix
        _coords (np.ndarray): Element coordinate matrix for graphical interface purposes
        gdl (np.ndarray): Degree of freedom matrix
    """

    def __init__(self, **kargs) -> None:
        """Create a continuum element

        Args:
            coords (np.ndarray): Element coordinate matrix
            _coords (np.ndarray): Element coordinate matrix for graphical interface purposes
            gdl (np.ndarray): Degree of freedom matrix
        """
        self.project_coords()
        self.project_disp()

    def project_coords(self) -> None:
        R = self.rotation_matrix()
        self.t_coords = self.coords@R
        self.t_coords = self.t_coords.reshape(
            len(self.t_coords), 1)  # Only works for bar elements

    def project_disp(self) -> None:
        """Project the displacements to the original coordinates system."""
        R = self.rotation_matrix()
        self.t_Ue = self.Ue.T@R
        # Only works for bar elements
        self.t_Ue = self.t_Ue.reshape(1, len(self.t_Ue))

    def setUe(self, U: np.ndarray) -> None:
        """Assing element local solution

        Args:
            U(np.ndarray): Global solution
        """

        for i in range(len(self.gdl)):
            self.Ue[i] = U[np.ix_(self.gdl[i])].flatten()
        n = len(self._coords)
        m = len(self.gdl)
        self._Ueg = self.Ue[np.ix_(np.linspace(
            0, m-1, m).astype(int), np.linspace(0, n-1, n).astype(int))]
        self._Ueg = np.array(self._Ueg.T.tolist()+[self._Ueg.T[0].tolist()])
        self.project_disp()

    def set_constitutive_model(self, CM) -> None:
        self.constitutive_model = CM

    def calculate_deformation_gradient(self, dpt):
        t_0 = self.t_coords.copy()
        t_t = t_0 + self.t_Ue.T
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

    def rotation_matrix(self, deformed=True) -> np.ndarray:
        """Calculate the rotation matrix for the element.

        This method should be implemented in derived classes.
        """
        raise NotImplementedError(
            "This method should be implemented in derived classes.")

    def transformation_matrix(self, deformed=True) -> np.ndarray:
        """Calculate the transformation matrix for the element.

        This method should be implemented in derived classes.
        """
        raise NotImplementedError(
            "This method should be implemented in derived classes.")

    def get_local_jacobian(self, jac: np.ndarray, deformed=True) -> np.ndarray:
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
        _j = self.jacs
        _dp = self.dpz
        weights = self.W

        Ke = 0.0
        Fe = 0.0
        for x, jac, wi, ni, dni in zip(_x, _j, weights, _p, _dp):
            J = self.get_local_jacobian(jac)
            detjac = np.linalg.det(J)
            # Esto asume que el jacobiano tiene inversa
            dpx = np.linalg.inv(J) @ dni
            # du = self.Ue.T @ dpx  # Creo
            F = self.calculate_deformation_gradient(dpx)
            E = self.green_lagrange_strain(F)
            C, S, const = self.constitutive_model(E)
            C = np.array([[C]])
            S_stiff, S_force = self.organize_S(S)

            BL = self.calculate_BL(dpx)
            BNL = self.calculate_BNL(dpx)
            Ke += const*(BL.T @ C @ BL + BNL.T @ S_stiff @ BNL) * detjac * wi
            Fe += const*(BL.T @ S_force) * detjac * wi

        T1 = self.transformation_matrix(False)
        T2 = self.transformation_matrix(False)
        Ke = T1.T @ Ke @ T1
        Fe = T2.T @ Fe
        return Ke, Fe
