if __name__ == '__main__':
    from FEM import NewtonTotalLagrangian, MGDCM, ContinumTotalLagrangian, QuadMembraneLinear, QuadShellLinear, Geometry3D
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.animation as animation
    coords = np.array([
        [0.0, 0.0, 0.0],   # Node 1
        [1.0, 0.0, 0.0],   # Node 2
        [1.0, 1.0, 0.1],   # Node 3 (slightly out of plane)
        [0.0, 1.0, 0.1]    # Node 4 (slightly out of plane)
    ])
    t = 0.1

    nu = 0.3
    young = 20500  # Young's modulus in MPa
    elements = [[0, 1, 2, 3]]
    types = [QuadShellLinear]*(len(elements))
    geo = Geometry3D(elements, coords, types, 5, fast=True)
    E = young  # Young's modulus in Pascals
    v = nu  # Poisson's ratio
    k = 5/6
    C = E/(1.0-v**2)*np.array([
        [1.0,   v, 0.0, 0.0, 0.0, 0.0],
        [v, 1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, (1.0-v)/2.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, (1.0-v)/2.0*k, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, (1.0-v)/2.0*k]])

    def cm(_E):
        E = np.zeros((6, 1))
        E[0, 0] = _E[0, 0]       # E11
        E[1, 0] = _E[1, 1]       # E22
        E[2, 0] = _E[2, 2]       # E33
        E[3, 0] = 2 * _E[0, 1]   # 2*E12
        E[4, 0] = 2 * _E[1, 2]   # 2*E23
        E[5, 0] = 2 * _E[0, 2]   # 2*E13

        S_voigt = C @ E
        S = np.zeros((3, 3))
        S[0, 0] = S_voigt[0, 0]  # S11
        S[1, 1] = S_voigt[1, 0]  # S22
        S[2, 2] = S_voigt[2, 0]  # S33
        S[0, 1] = S[1, 0] = S_voigt[3, 0]  # S12
        S[1, 2] = S[2, 1] = S_voigt[4, 0]  # S23
        S[0, 2] = S[2, 0] = S_voigt[5, 0]  # S13
        return C, S
    O = ContinumTotalLagrangian(
        geo, cm, solver=NewtonTotalLagrangian, override_nvn=True)
    np.random.seed(42)
    for e in O.elements:
        e3s = []
        for n in e.coords:
            vec = np.random.random(3)-0.5
            e3s.append(vec/np.linalg.norm(vec))
        e.set_thickness(t, e3s)
    e: QuadShellLinear = O.elements[0]
    W, _ = np.polynomial.legendre.leggauss(4)
    e.Ue = np.random.random(e.Ue.shape)*0.01  # small random displ
    deff = True
    for w in W:
        JS = e._get_spatial_derivatives(w, deformed=deff)
        for i, z in enumerate(e.Z):
            num = e.numerical_jacobian(*z, w, deformed=deff)
            delta = np.max(np.abs(JS[i]-num))
            if delta > 1e-6:
                print(f"Discrepancy at Gauss point {i}, w={w}: {delta}")
                print("Numerical Jacobian:")
                print(num)
                print("Analytical Jacobian:")
                print(JS)
            else:
                print(
                    f"Gauss point {i}, w={w}: Jacobian matches numerical calculation.")
    # New incremental displacements
    new_ue = np.random.random(e.Ue.shape)*0.1  # small random displ
    for w in W:
        JS, _JS = e.get_jacobians(w)
        FS = e.calculate_deformation_gradients(w)
        for i, z in enumerate(e.Z):
            JINV = _JS[i]
            dpx = e.calculate_dpxs(w, i, JINV)
            BNL = e.calculate_BNL(dpx)
            uk = new_ue.T.flatten()
            uij = BNL @ uk
            F = np.zeros([3, 3])
            F[0, 0] = uij[0] + 1
            F[1, 1] = uij[1] + 1
            F[2, 2] = uij[2] + 1
            F[0, 1] = uij[3]
            F[1, 0] = uij[4]
            F[0, 2] = uij[5]
            F[2, 0] = uij[6]
            F[1, 2] = uij[7]
            F[2, 1] = uij[8]
            UIJ = F - np.eye(3)
            uij_num = e.u_incr_deriv(*z, w, new_ue, h=1e-6)
            uij_num = JINV @ uij_num.T
            uij_num = uij_num.T
            a = 0
    for w in W:
        JS, _JS = e.get_jacobians(w)
        FS = e.calculate_deformation_gradients(w)
        thetas = e.derivatives_transformation()
        for i, z in enumerate(e.Z):
            JINV = _JS[i]
            dpx = e.calculate_dpxs(w, i, JINV)
            BNL = e.calculate_BNL(dpx)
            UIJ = FS[i] - np.eye(3)
            theta = thetas[i]
            a = 0
