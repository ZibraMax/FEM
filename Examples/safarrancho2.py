if __name__ == '__main__':
    from FEM import MGDCM, NewtonTotalLagrangian, ContinumTotalLagrangian, BarLinear, Geometry2D
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.animation import FuncAnimation

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
    h = 1
    t = 653
    A = h*t
    P = 1000
    young = 20500  # Young's modulus in MPa

    coords = [
        [0, 0],
        [a, b],
        [2*a, 0]
    ]
    elements = [[0, 1],
                [1, 2]]

    types = [BarLinear]*(len(elements))

    geo = Geometry2D(elements, coords, types, 2, fast=True)
    for node in [0, 2]:
        geo.cbe += [[node*2, 0], [node*2+1, 0]]

    # def cm(Ex):
    #     # Ogden hyperelastic constitutive model for bar elements
    #     C0 = young
    #     alfa = [3, 1]  # Linear

    #     pstr = np.real(np.sqrt(2 * Ex + 1))
    #     C0 = np.where(pstr < 1, 1 * C0, C0)

    #     Ct = C0 / (alfa[0] - alfa[1]) * ((alfa[0] - 2) * pstr **
    #                                      (alfa[0] - 4) - (alfa[1] - 2) * pstr**(alfa[1] - 4))
    #     Sx = C0 / (alfa[0] - alfa[1]) * \
    #         (pstr**(alfa[0] - 2) - pstr**(alfa[1] - 2))

    #     return Ct[0, 0], Sx[0, 0], A

    def cm(E):
        C = young
        S = young * E
        return C, S, A
    O = ContinumTotalLagrangian(geo, cm, solver=MGDCM, override_nvn=True)
    O.solver.set_delta_lambda_bar(0.1)
    O.solver.momentum = False
    O.solver.set_increments(220)
    O.cbn = [[3, -P]]
    O.solve()

    displacements = []
    load_factors = []
    for i in range(len(O.solver.solutions)):
        O.solver.setSolution(i, elements=True)
        displacements.append(-O.U[3][0])
        load_factors.append(O.solution_info['ld'])
    data = np.array([displacements, load_factors]).T
    np.savetxt('./Examples/examples_results/safarrancho2.csv', data,
               header='Displacement\tLoad factor', delimiter=',')
    us = np.linspace(0, np.max(displacements), len(displacements))
    force = analytical_force(us, young, A, a, b)

    # Plot the results
    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    cosa2, = ax2.plot(us[:1], -(force/P)[:1], '--', c='k',
                      lw=5, label="Analytical")
    cosa, = ax2.plot(
        displacements[:1], load_factors[:1], '-', label='Numerical')
    ax2.legend()

    ax2.set_ylabel('Load factor')
    ax2.set_xlabel('Displacement')
    ax2.set_title('Displacement vs Load factor')
    ax2.grid()
    ax2.set_xlim(1.1*min(displacements), 1.1*max(displacements))
    ax2.set_ylim(1.1*min(load_factors), 1.1*max(load_factors))
    ax.set(xlim=(-a*0.1, 2*a*1.1), xlabel='X')
    ax.set(ylim=(-b*1.5, b*1.1), ylabel='Y')
    plt.tight_layout()
    for e in O.elements:
        ax.plot(e.coords[:, 0], e.coords[:, 1], 'k-')
        L = np.linalg.norm(e.coords[1] - e.coords[0])
        print(L)
    # Deformed shape
    lines = []
    for e in O.elements:
        coords = e.coords + e.Ue.T
        lines.append(
            ax.plot(coords[:, 0], coords[:, 1], 'r-')[0])

    def animate(i, lines):
        O.solver.setSolution(i, elements=True)
        for j, e in enumerate(O.elements):
            coords = e.coords + e.Ue.T
            lines[j].set_data(coords.T)

        cosa.set_data(displacements[:i], load_factors[:i])
        cosa2.set_data(us[:i], -(force/P)[:i])
        return lines+[cosa, cosa2]

    anim = FuncAnimation(
        fig,
        animate,
        frames=len(O.solver.solutions)-1,
        interval=60,
        fargs=(lines,),
        blit=True
    )
    # html = anim.save('./Examples/examples_results/Truss_non_lineal.mp4')

    plt.show()
