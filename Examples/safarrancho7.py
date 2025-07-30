if __name__ == '__main__':
    from FEM import NewtonTotalLagrangian, MGDCM, ContinumTotalLagrangian, QuadMembraneLinear, QuadShellLinear, Geometry3D
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.animation as animation

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

    a = 400
    b = 20
    L = (a**2 + b**2)**0.5
    h = 653
    t = 1
    A = h*t
    P = 1000
    nu = 0.3
    young = 20500  # Young's modulus in MPa

    coords = [
        [0, 0, 0],
        [a, 0, b],
        [2*a, 0, 0],
        [0, h, 0],
        [a, h, b],
        [2*a, h, 0],
    ]

    elements = [
        [0, 1, 4, 3],
        [1, 2, 5, 4]]

    types = [QuadShellLinear]*(len(elements))

    geo = Geometry3D(elements, coords, types, 5, fast=True)

    disp_fixed_nodes = [0, 2, 3, 5]

    geo.cbe = []
    for i in disp_fixed_nodes:
        geo.cbe.append([i*5, 0])
        geo.cbe.append([i*5+1, 0])
        geo.cbe.append([i*5+2, 0])

    E = young  # Young's modulus in Pascals
    v = nu  # Poisson's ratio
    C = E/((1.0+v)*(1.0-2.0*v))*np.array([
        [1.0-v, v, v, 0.0, 0.0, 0.0],
        [v, 1.0-v, v, 0.0, 0.0, 0.0],
        [v, v, 1.0-v, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, (1.0-2.0*v)/2.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, (1.0-2.0*v)/2.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, (1.0-2.0*v)/2.0]])

    O = ContinumTotalLagrangian(
        geo, cm, solver=NewtonTotalLagrangian, override_nvn=True)
    for e in O.elements:
        e.set_thickness(t)
    O.solver.load_steps = 100

    nodes_force = [1, 4]

    O.cbn = [[nodes_force[0]*5+2, -P/2], [nodes_force[1]*5+2, -P/2]]
    O.solve()

    displacements = []
    load_factors = []
    for i in range(len(O.solver.solutions)):
        O.solver.setSolution(i, elements=True)
        displacements.append(-O.U[nodes_force[0]*5+2][0])
        load_factors.append(O.solution_info['ld'])
    data = np.array([displacements, load_factors]).T

    plots = []
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig.add_subplot(1, 2, 2)
    for e in O.elements:
        surf = ax.plot_trisurf(e.coords[:, 0], e.coords[:, 1],
                               e.coords[:, 2], alpha=0.5, color='b')
        plots.append(surf)
    pl, = ax2.plot(displacements[:1], load_factors[:1],
                   '-', c='k', label="Continumm incremental")
    plots.append(pl)
    ax2.set_xlabel('Displacement')
    ax2.set_ylabel('Load factor')

    def animate(i, plots):
        plots = []
        ax.clear()
        ax2.clear()
        O.solver.setSolution(i-1, elements=True)
        for e in O.elements:
            coords = e.coords + e.Ue[:3].T
            surf = ax.plot_trisurf(coords[:, 0], coords[:, 1],
                                   coords[:, 2], alpha=0.5, color='r')
            plots.append(surf)
        ax.set_xlim(0, 2*a)
        ax.set_ylim(0, h)
        ax.set_zlim(-b, b)
        ax.set_title(
            f'Deformed shape at load step {i}, load factor: {O.solution_info["ld"]:.2f}')
        plots.append(pl)
        pl, = ax2.plot(displacements[:i],
                       load_factors[:i], '-', c='k', label="Continumm incremental")
        plots.append(pl)
        ax2.legend()
        ax2.grid()
        ax2.set_xlabel('Displacement')
        ax2.set_ylabel('Load factor')
        return plots

    pam_ani = animation.FuncAnimation(fig, animate, fargs=(plots,),
                                      interval=5, blit=False, frames=len(O.solver.solutions))
    plt.show()
