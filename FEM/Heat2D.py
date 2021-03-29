"""2D Stady state heat with convective border conditions
"""


from .Core import Core, tqdm, np, Geometry
import matplotlib.pyplot as plt
from typing import Callable, Tuple


class Heat2D(Core):
    """Creates a Heat2D problem with convective borders

    The differential equation is:

    .. math::
        -\\frac{\\partial}{\\partial x}\\left(k_x\\frac{\\partial T}{\\partial x}\\right) - \\frac{\\partial}{\\partial y}\\left(k_y\\frac{\\partial T}{\\partial y}\\right)=f(x,y)

    With convective border conditions:

    .. math::
        k_x\\frac{\\partial T}{\\partial x}n_x+k_y\\frac{\\partial T}{\\partial y}n_y+\\beta (T-T_\\infty)=\\hat{q_n}

    Args:
        geometry (Geometry): Input 1 variable per node geometry
        kx (Tuple[float, list]): Heat transfer coeficient in x direction. If number, all element will have the same coefficient. If list, each position will be the element coefficient, so len(kx) == len(self.elements)
        ky (Tuple[float, list]): Heat transfer coeficient in y direction. If number, all element will have the same coefficient. If list, each position will be the element coefficient, so len(kx) == len(self.elements)
        Ta (float): Anbient temperature.
        f (Callable, optional): Internal heat generation function. Defaults to None.
        """

    def __init__(self, geometry: Geometry, kx: Tuple[float, list], ky: Tuple[float, list], Ta: float, f: Callable = None) -> None:
        """Creates a Heat2D problem with convective borders

        The differential equation is:

        .. math::
            -\\frac{\\partial}{\\partial x}\\left(k_x\\frac{\\partial T}{\\partial x}\\right) - \\frac{\\partial}{\\partial y}\\left(k_y\\frac{\\partial T}{\\partial y}\\right)=f(x,y)

        With convective border conditions:

        .. math::
            k_x\\frac{\\partial T}{\\partial x}n_x+k_y\\frac{\\partial T}{\\partial y}n_y+\\beta (T-T_\\infty)=\\hat{q_n}

        Args:
            geometry (Geometry): Input 1 variable per node geometry
            kx (Tuple[float, list]): Heat transfer coeficient in x direction. If number, all element will have the same coefficient. If list, each position will be the element coefficient, so len(kx) == len(self.elements)
            ky (Tuple[float, list]): Heat transfer coeficient in y direction. If number, all element will have the same coefficient. If list, each position will be the element coefficient, so len(kx) == len(self.elements)
            Ta (float): Anbient temperature.
            f (Callable, optional): Internal heat generation function. Defaults to None.
        """
        if type(kx) == float or type(kx) == int:
            kx = [kx]*len(geometry.elements)
        if type(ky) == float or type(ky) == int:
            ky = [ky]*len(geometry.elements)

        self.kx = kx
        self.ky = ky
        self.Ta = Ta
        self.f = f
        self.geometry = geometry
        Core.__init__(self, geometry)

    def elementMatrices(self) -> None:
        """Calculate the element matrices using Gauss Legendre quadrature.
        """
        for ee, e in enumerate(tqdm(self.elements, unit='Element')):
            m = len(e.gdl.T)
            K = np.zeros([m, m])
            H = np.zeros([m, m])
            F = np.zeros([m, 1])
            P = np.zeros([m, 1])
            _x, _p = e.T(e.Z.T)
            jac, dpz = e.J(e.Z.T)
            detjac = np.linalg.det(jac)
            _j = np.linalg.inv(jac)
            dpx = _j @ dpz

            for i in range(m):
                for j in range(m):
                    for k in range(len(e.Z)):
                        K[i, j] += (self.kx[ee]*dpx[k, 0, i]*dpx[k, 0, j] +
                                    self.ky[ee]*dpx[k, 1, i]*dpx[k, 1, j])*detjac[k]*e.W[k]
                if self.f:
                    for k in range(len(e.Z)):
                        F[i][0] += _p[k, i] * self.f(_x[k])*detjac[k]*e.W[k]

            if e.intBorders:
                for j in range(len(e.borders)):
                    border = e.borders[j]
                    if len(border.properties['load_x']):
                        _x, _p = e.T(e.Tj[j](border.Z.T))
                        _s = border.TS(border.Z.T)
                        detjac = border.coords[-1, 0]*0.5
                        for i in range(m):
                            for fx in border.properties['load_x']:
                                for k in range(len(border.Z)):
                                    P[i, 0] += self.Ta*fx(_s[k, 0])*_p[k, i] * \
                                        detjac*border.W[k]
                                for j in range(m):
                                    for k in range(len(border.Z)):
                                        H[i, j] += fx(_s[k, 0])*_p[k, i] * _p[k, j] * \
                                            detjac*border.W[k]
            e.Fe += F+P
            e.Ke += K+H

    def defineConvectiveBoderConditions(self, segment: int, beta: float = 0, tol: float = 1*10**(-5)) -> None:
        """Define convective borders

        Args:
            segment (int): Segment in wich load will be applied
            beta (float, optional): Convective coeficient :math:`\\beta` . Defaults to 0.
            tol (float, optional): Tolerancy to find adjacent nodes. Defaults to 1*10**(-5).
        """
        self.geometry.loadOnSegment(segment, fx=lambda s: beta, tol=tol)

    def postProcess(self) -> None:
        """Generate the temperature surface for the geometry

        """
        X = []
        Y = []
        U = []
        fig = plt.figure()
        ax = fig.add_subplot()
        for e in tqdm(self.elements, unit='Element'):
            _x, _u, _ = e.giveSolution(True)
            X += _x.T[0].tolist()
            Y += _x.T[1].tolist()
            U += _u[0].tolist()
        surf = ax.tricontourf(X, Y, U, cmap='magma', levels=20)
        plt.colorbar(surf, ax=ax)
        ax.set_title(r'$T$')
        mask = self.geometry.mask
        if self.geometry.holes:
            for hole in self.geometry.holes:
                Xs = np.array(hole['vertices'])[:, 0]
                Ys = np.array(hole['vertices'])[:, 1]
                ax.fill(Xs, Ys, color='white', zorder=30)
        if not mask == None:
            mask = np.array(mask)
            cornersnt = np.array(mask[::-1])
            xmin = np.min(cornersnt[:, 0])
            xmax = np.max(cornersnt[:, 0])
            ymin = np.min(cornersnt[:, 1])
            ymax = np.max(cornersnt[:, 1])
            Xs = [xmin, xmax, xmax, xmin]+cornersnt[:, 0].tolist()
            Ys = [ymin, ymin, ymax, ymax]+cornersnt[:, 1].tolist()
            ax.fill(Xs, Ys, color='white', zorder=30)
        ax.set_aspect('equal')
