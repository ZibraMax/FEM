import numpy as np


class Element():
    def __init__(self, coords, _coords, gdl):
        self.coords = coords
        self._coords = _coords
        self.gdl = gdl
        self.gdlm = []
        for i in range(len(self.gdl)):
            for j in range(len(self.gdl[i])):
                self.gdlm.append(self.gdl[i, j])
        self.n = int(len(self.gdl)*len(self.gdl[0]))
        self.Ke = np.zeros([self.n, self.n])
        self.Fe = np.zeros([self.n, 1])
        self.Ue = np.zeros(self.gdl.shape)
        self.Qe = np.zeros([self.n, 1])

    def T(self, z):
        p = self.psis(z)
        return p@self.coords, p

    def inverseMapping(self, x0, n=100):
        tol = 1*10**(-6)
        zi = np.zeros(x0.shape)+0.1
        for _ in range(n):
            xi = x0 - self.T(zi)[0].T
            _J = np.linalg.inv(self.J(zi)[0])
            xi = xi.T
            xi = xi.reshape(list(xi.shape)+[1])
            dz = _J@xi
            zi += dz[:, :, 0].T
            if np.max(np.abs(dz)) < tol:
                break
        return zi

    def J(self, z):
        dpsis = self.dpsis(z).T
        return dpsis @ self.coords, dpsis

    def giveSolution(self, SVSolution=False):
        _z = self.domain
        _x, _p = self.T(_z.T)
        if SVSolution:
            j, dpz = self.J(_z.T)  # TODO Revisar con Reddy
            dpx = np.linalg.inv(j) @ dpz
            # print((self.Ue @ np.transpose(dpx,axes=[0,2,1])).shape)
            # TODO REVISAR VS
            return _x, self.Ue@_p.T, self.Ue @ np.transpose(dpx, axes=[0, 2, 1])
        return _x, self.Ue@_p.T

    def giveSolutionPoint(self, Z, SVSolution=False):
        _x, _p = self.T(Z)
        if SVSolution:
            j, dpz = self.J(Z)  # TODO Revisar con Reddy
            dpx = np.linalg.inv(j) @ dpz
            # TODO REVISAR VS
            return _x, self.Ue@_p.T, self.Ue @ np.transpose(dpx, axes=[0, 2, 1])
        return _x, self.Ue@_p.T

    def setUe(self, U):
        for i in range(len(self.gdl)):
            self.Ue[i] = U[np.ix_(self.gdl[i])].flatten()
        n = len(self._coords)
        m = len(self.gdl)
        self._Ueg = self.Ue[np.ix_(np.linspace(
            0, m-1, m).astype(int), np.linspace(0, n-1, n).astype(int))]
        self._Ueg = np.array(self._Ueg.T.tolist()+[self._Ueg.T[0].tolist()])

    def integrate(self, f):
        integral = 0
        for w, z in zip(self.W, self.Z):
            integral += f(z)*w
        return integral
