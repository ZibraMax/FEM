from .LinealElement import LinealElement
from .LinearScheme import LinearScheme
import numpy as np
from typing import Callable


class OriHinge(LinealElement, LinearScheme):

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, **kargs):
        bar_coords = np.array([coords[1], coords[2]])
        print(bar_coords)
        LinealElement.__init__(self, bar_coords, gdl, **kargs)
        LinearScheme.__init__(self, 3)

        self.hinge_coords = coords
        self.calculate_vectors()
        self.theta_0 = self.calculate_theta()
        self.theta = self.calculate_theta()

    def calculate_vectors(self):
        self.i = self.hinge_coords[0] + self.Ue.T[0]
        self.j = self.hinge_coords[1] + self.Ue.T[1]
        self.k = self.hinge_coords[2] + self.Ue.T[2]
        self.l = self.hinge_coords[3] + self.Ue.T[3]

        self.rij = -(self.i-self.j)
        self.rkj = -(self.k-self.j)
        self.rkl = -(self.k-self.l)

        self.rij_n = np.linalg.norm(self.rij)
        self.rkj_n = np.linalg.norm(self.rkj)
        self.rkl_n = np.linalg.norm(self.rkl)

        self.m = np.cross(self.rij, self.rkj)
        self.n = np.cross(self.rkj, self.rkl)
        self.m_n2 = np.sum(self.m**2)
        self.n_n2 = np.sum(self.n**2)

    def calculate_theta(self):
        eta = 1
        if np.dot(self.m, self.rkl) != 0:
            eta = np.sign(np.dot(self.m, self.rkl))

        theta_ = eta*np.arccos(np.dot(self.m, self.n) /
                                     (np.linalg.norm(self.m)*np.linalg.norm(self.n)))
        theta_ = theta_ % (2*np.pi)
        return theta_

    def set_kf(self, kf: Callable):
        self.kf = kf

    def get_kf(self):
        return self.kf(self.theta)

    def jacobian(self):
        self.dtdxi = (self.rkj_n/self.m_n2) * self.m
        self.dtdxl = (-self.rkj_n/self.n_n2) * self.n
        self.dtdxj = (np.dot(self.rij, self.rkj)/self.rkj_n**2-1)*self.dtdxi - \
            np.dot(self.rkl, self.rkj)/self.rkj_n**2 * self.dtdxl

        self.dtdxk = (np.dot(self.rkl, self.rkj)/self.rkj_n**2-1)*self.dtdxl - \
            np.dot(self.rij, self.rkj)/self.rkj_n**2 * self.dtdxi
        J = np.array([self.dtdxi, self.dtdxj, self.dtdxk, self.dtdxl])
        return J

    def elementMatrix(self):
        self.calculate_vectors()
        self.theta = self.calculate_theta()
        J = self.jacobian().flatten().reshape([12, 1])
        k = self.get_kf()
        return k*J@J.T

    def elementMatrixNonLineal(self):
        J = self.jacobian().flatten().reshape([12, 1])
        k = self.get_kf()
        Ke = k*J@J.T
        Te = 0.0
        return Ke, Te
