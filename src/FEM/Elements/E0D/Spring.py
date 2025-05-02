import numpy as np


class Spring():
    def __init__(self, dof1, dof2, k):
        self.k = k
        self.dof1 = dof1
        self.dof2 = dof2
        self.gdl = np.array([[dof1, dof2]])
        self.internal_force = None
        self.theta0 = 210*np.pi/180  # Initial angle in radians

    def restartMatrix(self):
        self.Ue = np.zeros((2, 1))
        self.dof1 = int(self.dof1)
        self.dof2 = int(self.dof2)

    def setUe(self, Ue: np.ndarray):
        self.Ue = Ue.flatten()[[self.dof1, self.dof2]]

    def stiffness_matrix(self):
        return np.array([[self.k, -self.k],
                         [-self.k, self.k]])

    def force_vector(self):
        u1, u2 = self.Ue.flatten()
        deltau = self.theta0 + (self.Ue[1] - self.Ue[0])
        return np.array([self.k * (deltau), -self.k * (deltau)])

    def calculate_vectors(self):
        self.deltau = self.Ue[1] - self.Ue[0]

    def calculate_theta(self) -> float:
        return self.theta0 + self.deltau

    def elementMatrixNonLineal(self):
        K = self.stiffness_matrix()
        F = self.force_vector()
        return K, -F.reshape(2, 1)
