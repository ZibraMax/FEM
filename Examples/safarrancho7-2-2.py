if __name__ == '__main__':
    from FEM import NewtonTotalLagrangian, MGDCM, ContinumTotalLagrangian, QuadMembraneLinear, Geometry3D
    from FEM.Elements.ContinumElements import SerendipityShellLinear as SerendipityShell
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.animation as animation

    x1 = np.array([-1.0, -1.0,  0.0])*5
    x2 = np.array([1.0, -1.0,  0.0])*5
    x3 = np.array([1.0,  1.0,  0.0])*5
    x4 = np.array([-1.0,  1.0,  0.0])*5
    t = 0.5
    P = 7000
    nu = 0.3
    young = 20500  # Young's modulus in MPa

    coords = [x1, x2, x3, x4]
    # Create serendipity coordinates
    coords.append((x2+x1)/2)
    coords.append((x3+x2)/2)
    coords.append((x4+x3)/2)
    coords.append((x1+x4)/2)
    coords = np.array(coords)
    elements = [[0, 1, 2, 3, 4, 5, 6, 7]]

    types = [SerendipityShell]*(len(elements))

    geo = Geometry3D(elements, coords, types, 5, fast=True)

    disp_fixed_nodes = [0, 1]

    geo.cbe = []
    for i in disp_fixed_nodes:
        geo.cbe.append([i*5, 0])
        geo.cbe.append([i*5+1, 0])
        geo.cbe.append([i*5+2, 0])
        geo.cbe.append([i*5+3, 0])
        geo.cbe.append([i*5+4, 0])

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
        geo, cm, solver=MGDCM, override_nvn=True)
    for e in O.elements:
        e.set_thickness(t, n_gauss_thickness=1)
    O.solver.set_delta_lambda_bar(0.05)
    O.solver.momentum = False
    O.solver.set_increments(220)
    nodes_force = [2, 3, 6]
    O.cbn = [[i*5+1, P/len(nodes_force)] for i in nodes_force]
    O.solve()

    displacements = []
    load_factors = []
    for i in range(len(O.solver.solutions)):
        O.solver.setSolution(i, elements=True)
        displacements.append(O.U[nodes_force[0]*5+1][0])
        load_factors.append(O.solution_info['ld']*P/len(nodes_force))
    data = np.array([displacements, load_factors]).T
    np.savetxt("safarrppp2.csv", data, delimiter=",")
    plt.plot(*data.T)
    plt.grid()
    plt.show()
