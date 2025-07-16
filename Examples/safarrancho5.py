if __name__ == '__main__':
    from FEM import MGDCM, ContinumTotalLagrangian, TriMembraneLinear, BarLinear, Geometry3D
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.animation import FuncAnimation
    h1 = 1
    h2 = 0.5
    h = h1 + h2
    d1 = 2.5
    d2 = 4
    P = 5
    t = 0.1
    N = 200
    L0 = 1

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

    coords = [[0.0, 0.0, h]]
    gamma = 60 * np.pi / 180
    for i in range(0, 6):
        theta = i * gamma
        coords.append([d1*np.cos(theta), d1*np.sin(theta), h1])
    for i in range(0, 6):
        theta = i * gamma - np.pi/2
        coords.append([d2*np.cos(theta), d2*np.sin(theta), 0])
    coords = np.array(coords)

    elements = [[0, 1, 2],
                [0, 2, 3],
                [0, 3, 4],
                [0, 4, 5],
                [0, 5, 6],

                [5, 7, 6],
                [6, 8, 1],
                [1, 9, 2],
                [2, 10, 3],
                [3, 11, 4],
                [4, 12, 5],
                ]
    types = [TriMembraneLinear]*(len(elements))

    geo = Geometry3D(elements, coords, types, 3, fast=True)
    geo.cbe = []
    for node in [7, 8, 9, 10, 11, 12]:
        geo.cbe += [[node*3, 0], [node*3+1, 0], [node*3+2, 0]]
    for node in [0, 1, 2, 3, 4, 5, 6]:
        geo.cbe += [[node*3, 0], [node*3+1, 0]]
    geo.cbn = [[2, -P]]

    O = ContinumTotalLagrangian(geo, cm, solver=MGDCM, verbose=True)
    O.solver.set_delta_lambda_bar(L0)
    O.solver.momentum = False
    O.solver.tol = 1e-3
    O.solver.set_increments(N)
    O.solve()

    displacements = []
    load_factors = []
    for i in range(len(O.solver.solutions)):
        O.solver.setSolution(i, elements=True)
        displacements.append(-O.U[2][0])
        load_factors.append(O.solution_info['ld'])

    # Plot the results
    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig.add_subplot(1, 2, 2)
    cosa, = ax2.plot(
        displacements[:1], load_factors[:1], '-', label='Numerical')
    ax2.legend()
    ax2.set_ylabel('Load factor')
    ax2.set_xlabel('Displacement')
    ax2.grid()
    ax2.set_xlim(1.1*min(displacements), 4)
    ax2.set_ylim(-80, 70)
    ax.set(xlim3d=(-d1*1.2, d1*1.2), ylim3d=(-d1*1.2, d1*1.2),
           zlim3d=(0, h), xlabel='X', ylabel='Y', zlabel='Z')
    ax.set_aspect('equal')
    lines = []
    for e in O.elements:
        coords = e.coords + e.Ue.T
        lines.append(
            ax.plot(coords[:, 0], coords[:, 1], 'r-')[0])

    def animate(i, lines, cosa):
        O.solver.setSolution(i, elements=True)
        for j, e in enumerate(O.elements):
            coords = e.coords + e.Ue.T
            lines[j].set_data_3d(coords.T)
        cosa.set_data(displacements[:i], load_factors[:i])
        return lines + [cosa]

    anim = FuncAnimation(
        fig,
        animate,
        frames=len(O.solver.solutions)-1,
        interval=60,
        fargs=(lines, cosa),
        blit=True
    )
    # anim.save('./Examples/examples_results/RToff_non_lineal.mp4')

    plt.show()
