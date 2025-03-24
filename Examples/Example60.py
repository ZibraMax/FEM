if __name__ == '__main__':
    from FEM import BarAndHingeNonLinear
    from FEM.Geometry import Geometry3D
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib
    import json
    import logging
    matplotlib.rcParams['animation.embed_limit'] = 2**128
    from matplotlib.animation import FuncAnimation

    def R(theta):
        return theta*np.pi/180

    def D(theta):
        return theta*180/np.pi
    stifness_f = 10
    k_hinges = 1
    k_panels = stifness_f*k_hinges
    E = 1000
    A = 1

    # Miura ori coordinates
    miura_coords = np.loadtxt('./miura_coords.txt')
    # load the perimeter
    perimeter = np.loadtxt('./perimeter.txt').astype(int)-1
    # load the hinges
    hinges = np.loadtxt('./folds.txt').astype(int)-1
    # load the panels
    panels = np.loadtxt('./bending.txt').astype(int)-1

    hinges = hinges[:, [2, 0, 1, 3]]
    panels = panels[:, [2, 0, 1, 3]]

    hinge_bars = hinges[:, [1, 2]]
    panel_bars = panels[:, [1, 2]]

    elements = hinges.tolist() + panels.tolist() + perimeter.tolist() + \
        hinge_bars.tolist() + panel_bars.tolist()
    types = ['OH']*len(hinges) + ['OH']*len(panels) + ['L1V'] * \
        len(perimeter) + ['L1V']*len(hinge_bars) + ['L1V']*len(panel_bars)
    geo = Geometry3D(elements, miura_coords, types, 3, fast=True)
    algo = [[1, 1, 1, 1],
            [3, 1, 1, 1],
            [5, 1, 1, 1],
            [7, 1, 1, 1],
            [9, 1, 1, 1],
            [11, 1, 1, 1],
            [243, 1, 1, 1],
            [245, 1, 1, 1],
            [247, 1, 1, 1],
            [249, 1, 1, 1],
            [251, 1, 1, 1],
            [253, 1, 1, 1]]
    for e in algo:
        node = e[0]-1
        geo.cbe.append([node*3, 0.0])
        geo.cbe.append([node*3+1, 0.0])
        geo.cbe.append([node*3+2, 0.0])
    O = BarAndHingeNonLinear(geo, E, A, verbose=True)
    O.addLoadNode(126, [0.0, 0.0, -1.0])
    O.solver.set_increments(1000)
    O.solver.maxiter = 100
    O.solver.tol = 1e-3
    O.solver.set_delta_lambda_bar(0.5)

    for h in range(len(hinges)):
        O.elements[h].set_kf(k_hinges)
    for p in range(len(hinges), len(hinges)+len(panels)):
        O.elements[p].set_kf(k_panels)
    O.solve()
    O.exportJSON('./Examples/Mesh_tests/Bar_and_hinge_non_linear_loco.json')

    displacements = []
    load_factors = []
    for i in range(len(O.solver.solutions)):
        O.solver.setSolution(i, elements=True)
        displacements.append([-O.U[126*3+2][0]])
        load_factors.append([O.solution_info['ld']])
        print(
            f"Displacement: {-O.U[126*3+2][0]} Load factor: {O.solution_info['ld']}")
    displacements = np.array(displacements)
    load_factors = np.array(load_factors)
    # Plot the results
    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig.add_subplot(1, 2, 2)
    cosas = ax2.plot(displacements, load_factors, 'o-',
                     label="Disp")
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
    elementes = {}
    k = 0
    for i, e in enumerate(O.elements):
        if e.__class__.__name__ != 'OriHinge':
            coords = e.coords + e.Ue.T
            lines.append(
                ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], 'r-')[0])
            elementes[k] = i
            k += 1

    def animate(i, lines):
        O.solver.setSolution(i, elements=True)
        for j, k in elementes.items():
            e = O.elements[k]
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
    ax.axis('off')
    ax.set(xlim3d=(0, 30), xlabel='X')
    ax.set(ylim3d=(0, 30), ylabel='Y')
    ax.set(zlim3d=(0, 30), zlabel='Z')

    html = anim.to_jshtml()
    with open('./mori.html', 'w') as f:

        f.write(html)
