if __name__ == '__main__':
    from FEM import BarAndHingeNonLinear
    from FEM.Geometry import Geometry3D
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib
    matplotlib.rcParams['animation.embed_limit'] = 2**128
    from matplotlib.animation import FuncAnimation

    stifness_f = 10

    def R(theta):
        return theta*np.pi/180

    def D(theta):
        return theta*180/np.pi

    def k_hinges(theta):
        return 0.06

    def k_panels(theta):
        return stifness_f*k_hinges(theta)

    E = 1e4
    A = 0.2

    # Miura ori coordinates
    miura_coords = np.loadtxt('./miura_cell_coords.txt')
    elements = [[0, 1],  # Bar
                [1, 2],  # Bar
                [3, 4],  # Bar
                [4, 5],  # Bar
                [6, 7],  # Bar
                [7, 8],  # Bar
                [0, 3],  # Bar
                [3, 6],  # Bar
                [2, 5],  # Bar
                [5, 8],  # Bar
                [1, 4],  # Bar
                [4, 7],  # Bar
                [3, 1],  # Extra bars
                [1, 5],  # Extra bars
                [6, 4],  # Extra bars
                [4, 8],  # Extra bars
                [6, 4, 3, 1],  # Hinges
                [8, 4, 5, 1],  # Hinges
                [3, 1, 4, 5],  # Hinges
                [6, 4, 7, 8],  # Hinges
                [4, 3, 1, 0],  # Panels
                [4, 5, 1, 2],  # Panels
                [7, 6, 4, 3],  # Panels
                [7, 8, 4, 5],  # Panels
                ]

    types = ['L1V']*16 + ['OH']*8

    geo = Geometry3D(elements, miura_coords, types, 3, fast=True)
    ebc_x = [0, 3, 6]
    ebc_y = [0, 2]
    ebc_z = [0, 1, 2, 6, 7, 8]
    for node in ebc_x:
        geo.cbe += [[node*3, 0]]
    for node in ebc_y:
        geo.cbe += [[node*3+1, 0]]
    for node in ebc_z:
        geo.cbe += [[node*3+2, 0]]
    O = BarAndHingeNonLinear(geo, E, A, verbose=False)
    O.addLoadNode(2, [-1.0, 0.0, 0.0])
    O.addLoadNode(5, [-1.0, 0.0, 0.0])
    O.addLoadNode(8, [-1.0, 0.0, 0.0])
    O.solver.set_increments(1000)
    O.solver.maxiter = 20
    O.solver.set_delta_lambda_bar(0.1)
    hinges = [16, 17, 18, 19]
    panels = [20, 21, 22, 23]
    for h in hinges:
        O.elements[h].set_kf(k_hinges)
    for p in panels:
        O.elements[p].set_kf(k_panels)

    O.solve()
    O.exportJSON('./Examples/Mesh_tests/Bar_and_hinge_non_linear.json')

    displacements = []
    load_factors = []
    for i in range(len(O.solver.solutions)):
        O.solver.setSolution(i, elements=True)
        angles = []
        for hinge in hinges:
            O.elements[hinge].calculate_vectors()
            tt = D(O.elements[17].calculate_theta())
            angles.append(tt)
        print(O.solution_info['ld'])
        displacements.append(angles)
        load_factors.append([O.solution_info['ld']]*len(angles))
    displacements = np.array(displacements)
    load_factors = np.array(load_factors)
    # Plot the results
    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig.add_subplot(1, 2, 2)
    cosas = ax2.plot(displacements, load_factors, 'o-',
                     label=[str(i) for i in hinges])
    ax2.set_ylabel('Load factor')
    ax2.set_xlabel('Dihedral angle (degrees)')
    ax2.legend()
    ax2.grid()
    ax2.set_xlim(np.min(displacements), np.max(displacements))
    ax2.set_ylim(np.min(load_factors), np.max(load_factors))
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
        for j, cosa in enumerate(cosas):
            cosa.set_data(displacements[:i, j], load_factors[:i, j])
        return lines+cosas

    anim = FuncAnimation(
        fig,
        animate,
        frames=len(O.solver.solutions),
        interval=60,
        fargs=(lines,),
        blit=True
    )
    html = anim.to_jshtml()
    with open('./mori.html', 'w') as f:

        f.write(html)
