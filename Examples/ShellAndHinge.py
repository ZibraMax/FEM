if __name__ == '__main__':
    from FEM import ShellAndHingeNonLinear
    from FEM.Geometry import Geometry3D
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.animation import FuncAnimation

    def R(theta):
        return theta*np.pi/180

    def D(theta):
        return theta*180/np.pi

    kf = 0.0003
    L = 1
    EA = 10000*0.1

    x_coord = np.sin(np.pi/3)*L
    theta_0 = R(210)
    phi = 2*np.pi-theta_0
    z_coord = np.sin(phi)*x_coord
    x_coord2 = np.cos(phi)*x_coord

    coords = [[0, -0.5*L, 0],
              [0, 0.5*L, 0],
              [x_coord, 0, 0],
              [x_coord2, 0, z_coord],
              [0, 0, 0.0],
              [0.5*x_coord, 0.25*L, 0],
              [0.5*x_coord, -0.25*L, 0],
              [0.5*x_coord2, 0.25*L, 0.5*z_coord],
              [0.5*x_coord2, -0.25*L, 0.5*z_coord]]

    elements = [[0, 1, 2, 4, 5, 6],
                [1, 0, 3, 4, 8, 7],
                [2, 0, 1, 3]]

    types = ['MITC6']*(len(elements)-1) + ['OH']*1

    geo = Geometry3D(elements, coords, types, 5, fast=True)
    for node in [0, 1, 2, 4, 5, 6]:
        geo.cbe += [[node*5, 0], [node*5+1, 0], [node*5+2, 0]]
    O = ShellAndHingeNonLinear(geo, EA, 0.3, 1, verbose=False)
    O.solver.set_increments(50)
    O.solver.maxiter = 1000
    O.solver.tol = 1e-4
    O.solver.set_delta_lambda_bar(0.01)
    O.solver.max_iter_momentum = 15
    O.solver.min_iter_momentum = 1
    for e in O.elements:
        if e.__class__.__name__ == 'OriHinge':
            e.set_kf(kf)
    O.addLoadNode(3, [0.0, 0, 100])
    O.solve()
    O.exportJSON('./Examples/Mesh_tests/Shell_and_hinge_non_linear.json')

    displacements = []
    load_factors = []
    for i in range(len(O.solver.solutions)):
        O.solver.setSolution(i, elements=True)
        O.elements[-1].calculate_vectors()
        tt = D(2*np.pi-O.elements[-1].calculate_theta())
        print(tt, O.solution_info['ld'])
        displacements.append(tt)
        load_factors.append(O.solution_info['ld'])

    # Plot the results
    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig.add_subplot(1, 2, 2)
    cosa, = ax2.plot(
        displacements[:1], load_factors[:1], 'o-', label='Shell and hinge')
    data = np.loadtxt('./bah.csv', delimiter=',')
    displacements_bah = data[:, 1]
    load_factors_bah = data[:, 0]
    ax2.plot(displacements_bah, load_factors_bah, 'o-',
             label='Bar and hinge', color='orange')
    ax2.legend()
    ax2.set_ylabel('Load factor')
    ax2.set_xlabel('Hinge angle (degrees)')
    ax2.grid()
    ax2.set_xlim(min(displacements), 1.1*max(displacements))
    ax2.set_ylim(1.1*min(load_factors), 1.1*max(load_factors))
    ax.set(xlim3d=(-L*0.6, L*0.6), xlabel='X')
    ax.set(ylim3d=(-L*0.6, L*0.6), ylabel='Y')
    ax.set(zlim3d=(-L*0.6, L*0.6), zlabel='Z')

    plt.tight_layout()
    for e in O.elements:
        if e.__class__.__name__ != 'OriHinge':
            ax.plot_trisurf(e.coords[:, 0], e.coords[:, 1],
                            e.coords[:, 2], alpha=0.5, color='b')
            L = np.linalg.norm(e.coords[1] - e.coords[0])
    # Deformed shape
    lines = []
    for e in O.elements:
        if e.__class__.__name__ != 'OriHinge':
            coords = e.coords + e.Ue[:3].T
            ax.plot_trisurf(coords[:, 0], coords[:, 1],
                            coords[:, 2], triangles=[[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]], alpha=0.5, color='r')
            lines.append(
                ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], 'o')[0])

    def animate(i, lines):
        O.solver.setSolution(i-1, elements=True)
        for j, e in enumerate(O.elements):
            if e.__class__.__name__ != 'OriHinge':
                coords = e.coords + e.Ue[:3].T
                lines[j].set_data_3d(coords.T)

        cosa.set_data(displacements[:i], load_factors[:i])
        return lines+[cosa]

    # anim = FuncAnimation(
    #     fig,
    #     animate,
    #     frames=len(O.solver.solutions),
    #     interval=60,
    #     fargs=(lines,),
    #     blit=True
    # )
    # html = anim.save('./Examples/examples_results/Truss_non_lineal.mp4')
    animate(len(O.solver.solutions), lines)
    np.savetxt('./Pyuthon2.csv',
               np.array([displacements, load_factors]).T, delimiter=',')
    plt.show()
