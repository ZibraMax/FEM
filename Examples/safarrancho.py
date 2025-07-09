if __name__ == '__main__':
    from FEM import MGDCM, NewtonTotalLagrangian, ContinumTotalLagrangian, BarLinear, Geometry3D
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.animation import FuncAnimation

    def R(theta):
        return theta*np.pi/180

    def D(theta):
        return theta*180/np.pi

    L = 10.5
    young = 20000000
    A = 0.05*0.05
    EA = young*A

    delta_x = np.sin(R(60))*np.cos(R(45))*L
    delta_y = -np.cos(R(60))*np.cos(R(45))*L

    delta_z = np.sin(R(45))*L

    coords = [[0, 0, delta_z],
              [0, np.cos(R(45))*L, 0],
              [-delta_x, delta_y, 0],
              [delta_x, delta_y, 0]]

    elements = [[0, 1],
                [0, 2],
                [0, 3]]

    types = [BarLinear]*(len(elements))

    geo = Geometry3D(elements, coords, types, 3, fast=True)
    for node in [1, 2, 3]:
        geo.cbe += [[node*3, 0], [node*3+1, 0], [node*3+2, 0]]

    def cm(E):
        C = young
        S = young * E
        return C, S, A
    O = ContinumTotalLagrangian(geo, cm, solver=MGDCM)
    O.solver.set_delta_lambda_bar(0.5)
    O.solver.momentum = False
    O.solver.set_increments(250)
    O.cbn = [[2, -0.1*EA]]
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
    cosa, = ax2.plot(displacements[:1], load_factors[:1], 'o-')
    ax2.set_ylabel('Load factor')
    ax2.set_xlabel('Displacement')
    ax2.set_title('Displacement vs Load factor')
    ax2.grid()
    ax2.set_xlim(1.1*min(displacements), 1.1*max(displacements))
    ax2.set_ylim(1.1*min(load_factors), 1.1*max(load_factors))
    ax.set(xlim3d=(-L*0.6, L*0.6), xlabel='X')
    ax.set(ylim3d=(-L*0.6, L*0.6), ylabel='Y')
    ax.set(zlim3d=(-L*0.6, L*0.6), zlabel='Z')
    plt.tight_layout()
    for e in O.elements:
        ax.plot(e.coords[:, 0], e.coords[:, 1],
                e.coords[:, 2], 'k-')
        L = np.linalg.norm(e.coords[1] - e.coords[0])
        print(L)
    # Deformed shape
    lines = []
    for e in O.elements:
        coords = e.coords + e.Ue.T
        lines.append(
            ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], 'r-')[0])

    def animate(i, lines):
        O.solver.setSolution(i, elements=True)
        for j, e in enumerate(O.elements):
            coords = e.coords + e.Ue.T
            lines[j].set_data_3d(coords.T)

        cosa.set_data(displacements[:i], load_factors[:i])
        return lines+[cosa]

    anim = FuncAnimation(
        fig,
        animate,
        frames=len(O.solver.solutions),
        interval=60,
        fargs=(lines,),
        blit=True
    )
    # html = anim.save('./Examples/examples_results/Truss_non_lineal.mp4')
    plt.show()
