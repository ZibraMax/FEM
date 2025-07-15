if __name__ == '__main__':
    from FEM import MGDCM, NewtonTotalLagrangian, ContinumTotalLagrangian, QuadMembraneLinear, Geometry3D
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.animation as animation

    def analytical_force(u, E, A, L, h):
        L0 = np.sqrt(L**2 + h**2)
        y = h - u
        l = np.sqrt(L**2 + y**2)
        strain = (L0-l) / L0
        strain = strain + 1/2*strain**2
        N = E * A * strain
        vertical_component = (y / l)
        return -2 * N * vertical_component

    a = 400
    b = 20
    L = (a**2 + b**2)**0.5
    h = 653
    t = 1
    A = h*t
    P = 1000
    nu = 0.5
    young = 20500  # Young's modulus in MPa

    c11 = young / (1 - nu**2)
    c12 = nu * c11
    c22 = c11
    c66 = young / (2 * (1 + nu))
    C = np.array([
        [c11, c12, 0.0],
        [c12, c22, 0.0],
        [0.0, 0.0, c66]])

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

    types = [QuadMembraneLinear]*(len(elements))

    geo = Geometry3D(elements, coords, types, 3, fast=True)
    geo.cbe = [[0, 0],
               [1, 0],
               [2, 0],
               [6, 0],
               [7, 0],
               [8, 0],
               [9, 0],
               [11, 0],
               [15, 0],
               [17, 0]]

    def cm(_E):
        E = np.zeros((3, 1))
        E[0, 0] = _E[0, 0]
        E[1, 0] = _E[1, 1]
        E[2, 0] = 2*_E[0, 1]
        S = C @ E
        _S = np.zeros((2, 2))
        _S[0, 0] = S[0, 0]
        _S[1, 1] = S[1, 0]
        _S[0, 1] = S[2, 0]
        _S[1, 0] = S[2, 0]
        return C, _S, t
    O = ContinumTotalLagrangian(geo, cm, solver=MGDCM, override_nvn=True)
    O.solver.set_delta_lambda_bar(0.05)
    O.solver.momentum = False
    O.solver.set_increments(220)
    O.cbn = [[5, -P/2], [14, -P/2]]
    O.solve()

    displacements = []
    load_factors = []
    for i in range(len(O.solver.solutions)):
        O.solver.setSolution(i, elements=True)
        displacements.append(-O.U[5][0])
        load_factors.append(O.solution_info['ld'])
    data = np.array([displacements, load_factors]).T
    np.savetxt('./Examples/examples_results/safarrancho3.csv', data,
               header='Displacement\tLoad factor', delimiter=',')
    us = np.linspace(0, np.max(displacements), len(displacements))
    force = -analytical_force(us, young, A, a, b)/P

    plots = []
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig.add_subplot(1, 2, 2)
    for e in O.elements:
        surf = ax.plot_trisurf(e.coords[:, 0], e.coords[:, 1],
                               e.coords[:, 2], alpha=0.5, color='b')
        plots.append(surf)
    pl, = ax2.plot(us[:1], force[:1], '--', lw=3, c='gray', label="Analytical")
    plots.append(pl)
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
        pl, = ax2.plot(us[:i], force[:i], '--', c='gray', lw=3,
                       label="Analytical")
        plots.append(pl)
        pl, = ax2.plot(displacements[:i],
                       load_factors[:i], '-', c='k', label="Continumm incremental")
        plots.append(pl)
        ax2.legend()
        ax2.grid()
        ax2.set_xlabel('Displacement')
        ax2.set_ylabel('Load factor')
        ax2.set_xlim(0, 60)
        ax2.set_ylim(-1, 1)
        return plots

    pam_ani = animation.FuncAnimation(fig, animate, fargs=(plots,),
                                      interval=5, blit=False, frames=len(O.solver.solutions))
    pam_ani.save(
        './Examples/examples_results/Truss_non_lineal_continumm_incremental_membranes.mp4')
    plt.show()
