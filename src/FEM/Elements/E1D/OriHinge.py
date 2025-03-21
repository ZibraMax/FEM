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
        self.A = self.rij@self.rkj/self.rkj_n**2
        self.B = self.rkl@self.rkj/self.rkj_n**2
        self.dA_dxj = (1 / self.rkj_n**2) * \
            ((2 * self.A - 1) * self.rkj - self.rij)
        self.dB_dxj = (1 / self.rkj_n**2) * (2 * self.B * self.rkj - self.rkl)
        self.dA_dxk = (1 / self.rkj_n**2) * (-2 * self.A * self.rkj + self.rij)
        self.dB_dxk = (1 / self.rkj_n**2) * \
            ((1 - 2 * self.B) * self.rkj + self.rkl)

        self.dtdxi = (self.rkj_n/self.m_n2) * self.m
        self.dtdxl = (-self.rkj_n/self.n_n2) * self.n
        self.dtdxj = (np.dot(self.rij, self.rkj)/self.rkj_n**2-1)*self.dtdxi - \
            np.dot(self.rkl, self.rkj)/self.rkj_n**2 * self.dtdxl

        self.dtdxk = (np.dot(self.rkl, self.rkj)/self.rkj_n**2-1)*self.dtdxl - \
            np.dot(self.rij, self.rkj)/self.rkj_n**2 * self.dtdxi

        cross_product = np.cross(self.rkj, self.m)
        d2tdxi2 = -(self.rkj_n / (self.m_n2**2)) * \
            self.dia(self.m, cross_product)

        cross_product_n = np.cross(self.rkj, self.n)
        d2tdxl2 = (self.rkj_n / (self.n_n2**2)) * \
            self.dia(self.n, cross_product_n)

        cross_product_m_rij = np.cross(self.rij, self.m)
        d2tdxixk = (np.outer(self.m, self.rkj) / (self.m_n2 * self.rkj_n)) + \
            (self.rkj_n / (self.m_n2**2)) * \
            self.dia(self.m, cross_product_m_rij)

        cross_product_n_rkl = np.cross(self.rkl, self.n)
        d2tdxlxj = (np.outer(self.n, self.rkj) / (self.n_n2 * self.rkj_n)) - \
            (self.rkj_n / (self.n_n2**2)) * \
            self.dia(self.n, cross_product_n_rkl)

        cross_product_m_diff = np.cross(self.rkj - self.rij, self.m)
        d2tdxixj = -(np.outer(self.m, self.rkj) / (self.m_n2 * self.rkj_n)) + \
            (self.rkj_n / (self.m_n2**2)) * \
            self.dia(self.m, cross_product_m_diff)

        cross_product_n_diff = np.cross(self.rkj - self.rkl, self.n)
        d2tdxlxk = -(np.outer(self.n, self.rkj) / (self.n_n2 * self.rkj_n)) - \
            (self.rkj_n / (self.n_n2**2)) * \
            self.dia(self.n, cross_product_n_diff)

        d2tdxj2 = np.outer(self.dtdxi, self.dA_dxj) + (self.A - 1) * d2tdxixj - (
            np.outer(self.dtdxl, self.dB_dxj) + self.B * d2tdxlxk
        )

        d2tdxjxk = np.outer(self.dtdxi, self.dA_dxk) + (self.A - 1) * d2tdxixk - (
            np.outer(self.dtdxl, self.dB_dxk) + self.B * d2tdxlxk
        )

        d2tdxk2 = np.outer(self.dtdxl, self.dB_dxk) + (self.B - 1) * d2tdxlxk - (
            np.outer(self.dtdxi, self.dA_dxk) + self.A * d2tdxixk
        )

        d2tdxlxi = np.zeros((3, 3))

        # 12x12 hesian matrix

        self.d2theta_dxi2 = np.block([[d2tdxi2, d2tdxixj, d2tdxixk, d2tdxlxi],
                                      [d2tdxixj.T, d2tdxj2, d2tdxjxk, d2tdxlxj.T],
                                      [d2tdxixk.T, d2tdxjxk.T,
                                          d2tdxk2, d2tdxlxk.T],
                                      [d2tdxlxi.T, d2tdxlxj, d2tdxlxk, d2tdxl2]])

    def dia(self, a, b):
        return np.outer(a, b) + np.outer(b, a)

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

    def get_M(self):
        return self.kf(self.theta_0)*(self.theta-self.theta_0)

    def jacobian(self):
        J = np.array([self.dtdxi, self.dtdxj, self.dtdxk, self.dtdxl])
        return J

    def elementMatrix(self):
        J = self.jacobian().flatten().reshape([12, 1])
        k = self.get_kf()
        return k*J@J.T

    def elementMatrixNonLineal(self):
        self.calculate_vectors()
        self.theta = self.calculate_theta()
        J = self.jacobian().flatten().reshape([12, 1])
        k = self.get_kf()
        Ke = k*J@J.T
        M = self.get_M()
        Kg = M*self.d2theta_dxi2.T
        Te = M*J
        return Ke+Kg, Te

        # self.calculate_vectors()
        # self.theta = self.calculate_theta()
