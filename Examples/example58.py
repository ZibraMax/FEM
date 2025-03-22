if __name__ == '__main__':
    from FEM import BarAndHingeNonLinear
    from FEM.Geometry import Geometry3D
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.animation import FuncAnimation

    def R(theta):
        return theta*np.pi/180

    def D(theta):
        return theta*180/np.pi

    def kf(theta):
        return 1

    L = 1
    EA = 100

    x_coord = np.sin(np.pi/3)*L
    theta_0 = R(210)
    phi = 2*np.pi-theta_0
    z_coord = np.sin(phi)*x_coord
    x_coord2 = np.cos(phi)*x_coord

    coords = ([[0, -0.5*L, 0],
               [0, 0.5*L, 0],
               [x_coord, 0, 0],
               [x_coord2, 0, z_coord]])

    elements = [[0, 1],
                [1, 2],
                [2, 0],
                [0, 3],
                [1, 3],
                [3, 0, 1, 2]]

    types = ['L1V']*(len(elements)-1) + ['OH']*1

    geo = Geometry3D(elements, coords, types, 3, fast=True)
    for node in [0, 1, 2]:
        geo.cbe += [[node*3, 0], [node*3+1, 0], [node*3+2, 0]]
    O = BarAndHingeNonLinear(geo, EA, 1, verbose=False)
    O.solver.set_increments(100)
    O.solver.maxiter = 500
    O.solver.tol = 1e-3
    O.solver.set_delta_lambda_bar(0.005)
    for e in O.elements:
        if e.__class__.__name__ == 'OriHinge':
            e.set_kf(kf)
    O.addLoadNode(3, [0.0, 0, 1])
    O.solve()
    O.exportJSON('./Examples/Mesh_tests/Bar_and_hinge_non_linear.json')

    displacements = []
    load_factors = []
    for i in range(len(O.solver.solutions)):
        O.solver.setSolution(i, elements=True)
        O.elements[-1].calculate_vectors()
        tt = D(O.elements[-1].calculate_theta())
        print(O.solution_info['ld'])
        displacements.append(tt)
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
    ax2.set_xlim(min(displacements), 1.1*max(displacements))
    ax2.set_ylim(1.1*min(load_factors), 1.1*max(load_factors))
    ax.set(xlim3d=(-L*0.6, L*0.6), xlabel='X')
    ax.set(ylim3d=(-L*0.6, L*0.6), ylabel='Y')
    ax.set(zlim3d=(-L*0.6, L*0.6), zlabel='Z')
    plt.tight_layout()
    for e in O.elements:
        if e.__class__.__name__ != 'OriHinge':
            ax.plot(e.coords[:, 0], e.coords[:, 1],
                    e.coords[:, 2], 'k-')
            L = np.linalg.norm(e.coords[1] - e.coords[0])
    # Deformed shape
    lines = []
    for e in O.elements:
        if e.__class__.__name__ != 'OriHinge':
            coords = e.coords + e.Ue.T
            lines.append(
                ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], 'r-')[0])

    def animate(i, lines):
        O.solver.setSolution(i, elements=True)
        for j, e in enumerate(O.elements):
            if e.__class__.__name__ != 'OriHinge':
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
    html = anim.save('./Examples/examples_results/Truss_non_lineal.mp4')
    plt.show()
