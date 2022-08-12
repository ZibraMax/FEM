"""1D Stady state heat with convective border conditions
"""


from .Core import Core, tqdm, np, Geometry
import matplotlib.pyplot as plt


class Heat1D(Core):
    """Creates a 1D Stady state heat problem. Convective border conditions can be applied

    The differential equation is:

    .. math::
        -\\frac{d}{dx}\\left(Ak\\frac{dT}{dx}\\right)+\\beta P(T-T_{\\infty})=0

    Args:
        geometry (Geometry): Input 1D Geometry. 1 variable per node
        A (float or list): Section area. If float, all elements will have the same area. If list, the i-element will have the i-area
        P (float or list): Section perimeter. If float, all elements will have the same perimeter. If list, the i-element will have the i-perimeter
        k (float or list): Conductivity. If float, all elements will have the same conductivity. If list, the i-element will have the i-conductivity
        beta (float or list): Transfer coeficient. If float, all elements will have the same transfer coeficient. If list, the i-element will have the i-transfer coeficient
        Ta (float): Ambient temperature (also called T∞)
        q (float or list, optional): Internal heat generation rate. If float, all elements will have the same internal heat generation rate coeficient. If list, the i-element will have the i-internal heat generation rate. Defaults to 0.0.
    """

    def __init__(self, geometry: Geometry, A: float, P: float, ku: float, beta: float, Ta: float, q: float = 0.0) -> None:
        """Creates a 1D Stady state heat problem. Convective border conditions can be applied

        The differential equation is:

        .. math::
            -\\frac{d}{dx}\\left(Ak\\frac{dT}{dx}\\right)+\\beta P(T-T_{\\infty})

        Args:
            geometry (Geometry): Input 1D Geometry. 1 variable per node
            A (float or list): Section area. If float, all elements will have the same area. If list, the i-element will have the i-area
            P (float or list): Section perimeter. If float, all elements will have the same perimeter. If list, the i-element will have the i-perimeter
            ku (float or list): Conductivity. If float, all elements will have the same conductivity. If list, the i-element will have the i-conductivity
            beta (float or list): Transfer coeficient. If float, all elements will have the same transfer coeficient. If list, the i-element will have the i-transfer coeficient
            Ta (float): Ambient temperature
            q (float or list, optional): Internal heat generation rate. If float, all elements will have the same internal heat generation rate coeficient. If list, the i-element will have the i-internal heat generation rate. Defaults to 0.0.
        """
        if isinstance(A, float) or isinstance(A, int):
            A = [A]*len(geometry.elements)
        if isinstance(P, float) or isinstance(P, int):
            P = [P]*len(geometry.elements)
        if isinstance(ku, float) or isinstance(ku, int):
            ku = [ku]*len(geometry.elements)
        if isinstance(beta, float) or isinstance(beta, int):
            beta = [beta]*len(geometry.elements)
        if isinstance(q, float) or isinstance(q, int):
            q = [q]*len(geometry.elements)
        self.A = A
        self.P = P
        self.ku = ku
        self.beta = beta
        self.Ta = Ta
        self.q = q
        Core.__init__(self, geometry)
        self.name = '1D Heat transfer'
        self.properties['A'] = self.A
        self.properties['P'] = self.P
        self.properties['ku'] = self.ku
        self.properties['beta'] = self.beta
        self.properties['Ta'] = self.Ta
        self.properties['q'] = self.q

    def elementMatrices(self) -> None:
        """Calculate the element matrices using Gauss Legendre quadrature.
        """

        for ee, e in enumerate(tqdm(self.elements, unit='Element')):
            m = len(e.gdl.T)
            K = np.zeros([m, m])
            F = np.zeros([m, 1])
            _x, _p = e.T(e.Z.T)
            jac, dpz = e.J(e.Z.T)
            detjac = np.linalg.det(jac)
            _j = np.linalg.inv(jac)
            dpx = _j @ dpz
            for i in range(m):
                for j in range(m):
                    for k in range(len(e.Z)):
                        K[i, j] += (self.ku[ee]*self.A[ee]*dpx[k][0][i]*dpx[k][0]
                                    [j]+self.beta[ee]*self.P[ee]*_p[k][i]*_p[k][j])*detjac[k]*e.W[k]
                for k in range(len(e.Z)):
                    F[i][0] += (_p[k][i]*(self.A[ee]*self.q[ee] +
                                self.P[ee]*self.beta[ee]*self.Ta))*detjac[k]*e.W[k]
            e.Fe += F
            e.Ke += K

    def defineConvectiveBoderConditions(self, node: int, value: float = 0) -> None:
        """Add a convective border condition. The value is: :math:`kA\\frac{dT}{dx}+\\beta A(T-T_{\infty})=value`

        Args:
            node (int): Node where the above border condition is applied
            value (float, optional): Defined below. Defaults to 0.
        """

        near = np.infty
        for i, e in enumerate(self.elements):
            act = min(abs(self.geometry.gdls[node][0] - e._coords[0]),
                      abs(self.geometry.gdls[node][0] - e._coords[0]))
            if act < near:
                near = act
                k = i
            if act == 0:
                break
        self.cbn += [[node, value+self.Ta*self.beta[k]*self.A[k]]]
        self.K[node, node] += self.beta[k]*self.A[k]

    def postProcess(self) -> None:
        """Generate graph of solution and solution derivative
        """

        X = []
        U1 = []
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 1, 1)
        for e in tqdm(self.elements, unit='Element'):
            _x, _u = e.giveSolution(False)
            X += _x.T[0].tolist()
            U1 += _u[0].tolist()
        ax1.plot(X, U1)
        ax1.grid()
        ax1.set_title(r'$T(x)$')
        plt.show()
