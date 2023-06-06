"""1D Stady state heat with convective border conditions
"""


from .Core import Core, tqdm, np, Geometry, CoreParabolic
import matplotlib.pyplot as plt


class Heat1D(Core):
    """Creates a 1D Stady state heat problem. Convective border conditions can be applied

        The differential equation is:

        .. math::
            -\\frac{d}{dx}\\left(Ak\\frac{dT}{dx}\\right)+\\beta P(T-T_{\\infty})=q

        Args:
            geometry (Geometry): Input 1D Geometry. 1 variable per node
            A (float or list): Section area. If float, all elements will have the same area. If list, the i-element will have the i-area
            P (float or list): Section perimeter. If float, all elements will have the same perimeter. If list, the i-element will have the i-perimeter
            ku (float or list): Conductivity. If float, all elements will have the same conductivity. If list, the i-element will have the i-conductivity
            beta (float or list): Transfer coeficient. If float, all elements will have the same transfer coeficient. If list, the i-element will have the i-transfer coeficient
            Ta (float): Ambient temperature
            q (float or list, optional): Internal heat generation rate. If float, all elements will have the same internal heat generation rate coeficient. If list, the i-element will have the i-internal heat generation rate. Defaults to 0.0.
        """

    def __init__(self, geometry: Geometry, A: float, P: float, ku: float, beta: float, Ta: float, q: float = 0.0) -> None:
        """Creates a 1D Stady state heat problem. Convective border conditions can be applied

        The differential equation is:

        .. math::
            -\\frac{d}{dx}\\left(Ak\\frac{dT}{dx}\\right)+\\beta P(T-T_{\\infty})=q

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


class Heat1DTransient(CoreParabolic):
    """Creates a 1D Stady state heat problem. Convective border conditions can be applied

    The differential equation is:

    .. math::
        \\frac{\\partial T}{\\partial t}-\\frac{\\partial}{\\partial x}\\left(Ak\\frac{\\partial T}{\\partial x}\\right)+\\beta P(T-T_{\\infty})=q

    Args:
        geometry (Geometry): Input 1D Geometry. 1 variable per node
        A (float or list): Section area. If float, all elements will have the same area. If list, the i-element will have the i-area
        P (float or list): Section perimeter. If float, all elements will have the same perimeter. If list, the i-element will have the i-perimeter
        k (float or list): Conductivity. If float, all elements will have the same conductivity. If list, the i-element will have the i-conductivity
        beta (float or list): Transfer coeficient. If float, all elements will have the same transfer coeficient. If list, the i-element will have the i-transfer coeficient
        Ta (float): Ambient temperature (also called T∞)
        q (float or list, optional): Internal heat generation rate. If float, all elements will have the same internal heat generation rate coeficient. If list, the i-element will have the i-internal heat generation rate. Defaults to 0.0.
    """

    def __init__(self, geometry: Geometry, A: float, P: float, ku: float, beta: float, Ta: float, q: float = 0.0, **kargs):
        """Creates a 1D Stady state heat problem. Convective border conditions can be applied

        The differential equation is:

        .. math::
            \\frac{\\partial T}{\\partial t}-\\frac{\\partial}{\\partial x}\\left(Ak\\frac{\\partial T}{\\partial x}\\right)+\\beta P(T-T_{\\infty})=q

        Args:
            geometry (Geometry): Input 1D Geometry. 1 variable per node
            A (float or list): Section area. If float, all elements will have the same area. If list, the i-element will have the i-area
            P (float or list): Section perimeter. If float, all elements will have the same perimeter. If list, the i-element will have the i-perimeter
            k (float or list): Conductivity. If float, all elements will have the same conductivity. If list, the i-element will have the i-conductivity
            beta (float or list): Transfer coeficient. If float, all elements will have the same transfer coeficient. If list, the i-element will have the i-transfer coeficient
            Ta (float): Ambient temperature (also called T∞)
            q (float or list, optional): Internal heat generation rate. If float, all elements will have the same internal heat generation rate coeficient. If list, the i-element will have the i-internal heat generation rate. Defaults to 0.0.
        """
        CoreParabolic.__init__(
            self, geometry=geometry, **kargs)

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
        self.name = '1D Heat transfer'
        self.properties['A'] = self.A
        self.properties['P'] = self.P
        self.properties['ku'] = self.ku
        self.properties['beta'] = self.beta
        self.properties['Ta'] = self.Ta
        self.properties['q'] = self.q
        self.convective_conditions = []

    def elementMatrices(self) -> None:
        """Calculate the element matrices using Gauss Legendre quadrature.
        """
        a1 = self.alpha*self.dt
        a2 = (1-self.alpha)*self.dt

        for ee, e in enumerate(tqdm(self.elements, unit='Element')):
            m = len(e.gdl.T)
            M1 = np.zeros([m, m])
            Kt = np.zeros([m, m])
            Ft = np.zeros([m, 1])
            Ktp1 = np.zeros([m, m])
            Ftp1 = np.zeros([m, 1])
            _x, _p = e.T(e.Z.T)
            jac, dpz = e.J(e.Z.T)
            detjac = np.linalg.det(jac)
            _j = np.linalg.inv(jac)
            dpx = _j @ dpz
            for i in range(m):
                for j in range(m):
                    for k in range(len(e.Z)):
                        Kt[i, j] += (self.ku[ee]*self.A[ee]*dpx[k][0][i]*dpx[k][0]
                                     [j]+self.beta[ee]*self.P[ee]*_p[k][i]*_p[k][j])*detjac[k]*e.W[k]

                        Ktp1[i, j] += (self.ku[ee]*self.A[ee]*dpx[k][0][i]*dpx[k][0]
                                       [j]+self.beta[ee]*self.P[ee]*_p[k][i]*_p[k][j])*detjac[k]*e.W[k]

                        M1[i, j] += (_p[k][i]*_p[k][j])*detjac[k]*e.W[k]
                for k in range(len(e.Z)):
                    Ft[i][0] += (_p[k][i]*(self.A[ee]*self.q[ee] +
                                           self.P[ee]*self.beta[ee]*self.Ta))*detjac[k]*e.W[k]
                    Ftp1[i][0] += (_p[k][i]*(self.A[ee]*self.q[ee] +
                                             self.P[ee]*self.beta[ee]*self.Ta))*detjac[k]*e.W[k]
            Fttp1 = self.dt*(self.alpha*Ft + (1-self.alpha)*Ftp1)
            e.Fe = (M1 - a2*Kt)@e.Ue.T + Fttp1
            e.Ke = M1 + a1*Ktp1

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
        self.convective_conditions.append([node, self.beta[k]*self.A[k]])

    def set_convective_conditions(self):
        for idx, value in self.convective_conditions:
            self.K[idx, idx] += value

    def borderConditions(self) -> None:
        a = super().borderConditions()
        self.set_convective_conditions()
        return a

    def postProcess(self, node=None, t0=None, tf=None, steps=None, dt=None, ax=None) -> None:
        """Post process the solution and steps
        """
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
        steps = len(self.solver.solutions)
        if not node:
            for i in range(steps):
                self.solver.setSolution(i, True)
                X = self.geometry.gdls.flatten().tolist()
                U1 = self.U.flatten().tolist()
                color = "gray"
                if i == 0:
                    color = "black"
                elif i == steps:
                    color = "red"

                ax.plot(X, U1, color=color)
            ax.set_xlabel("X")
            ax.set_ylabel(f"T all nodes")
        else:
            color = "k"
            X = []
            U1 = []
            for i in range(steps):
                self.solver.setSolution(i, True)
                U1.append(self.U[node])
                X.append(self.solver.solutions_info[i]["time"])

            ax.plot(X, U1, color=color)
            ax.set_xlabel("t")
            ax.set_ylabel(f"T node {node}")
            return X, U1
        ax.grid()
        ax.set_title(r'$T(x)$')
