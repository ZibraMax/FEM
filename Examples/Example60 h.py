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

    data = {
        "vertices": [
            [0.0, 0.0, 0.0],
            [5.0, 2.0, 0.0],
            [10.0, 0.0, 0.0],
            [0.0, 5.0, 3.0],
            [5.0, 7.0, 3.0],
            [10.0, 5.0, 3.0],
            [0.0, 10.0, 0.0],
            [5.0, 12.0, 0.0],
            [10.0, 10.0, 0.0]
        ],
        "lines": [
            [0, 1],
            [1, 2],
            [3, 4],
            [4, 5],
            [6, 7],
            [7, 8],
            [1, 4],
            [4, 7],
            [0, 3],
            [3, 6],
            [2, 5],
            [5, 8],
            [1, 3],
            [1, 5],
            [4, 6],
            [4, 8]
        ],
        "bending": [
            [0, 3, 1, 4],
            [2, 1, 5, 4],
            [3, 6, 4, 7],
            [5, 4, 8, 7]
        ],
        "folds": [
            [3, 1, 4, 5],
            [6, 4, 7, 8],
            [1, 3, 4, 6],
            [1, 4, 5, 8]
        ],
        "ebc": [
            [0, 1, 1, 1],
            [2, 0, 1, 1],
            [3, 1, 0, 0],
            [6, 1, 0, 1],
            [8, 0, 0, 1]
        ],
        "nbc": [
            [2, -1.0, 0, 0],
            [8, -1.0, 0, 0]
        ]
    }

    hinges = data["folds"]
    panels = data["bending"]
    lines = data["lines"]
    miura_coords = np.array(data["vertices"])

    elements = data["folds"] + data["bending"] + data["lines"]

    types = ['OH']*len(hinges) + ['OH']*len(panels) + ['L1V']*len(lines)
    geo = Geometry3D(elements, miura_coords, types, 3, fast=True)
    algo = data["ebc"]
    for e in algo:
        node = e[0]
        if e[1] == 1:
            geo.cbe.append([node*3, 0.0])
        if e[2] == 1:
            geo.cbe.append([node*3+1, 0.0])
        if e[3] == 1:
            geo.cbe.append([node*3+2, 0.0])
    O = BarAndHingeNonLinear(geo, E, A, verbose=False)
    O.addLoadNode(2, [-1.0, 0.0, 0.0])
    O.addLoadNode(8, [-1.0, 0.0, 0.0])
    O.solver.set_increments(100)
    O.solver.maxiter = 100
    O.solver.tol = 1e-5
    O.solver.set_delta_lambda_bar(0.5)

    for h in range(len(hinges)):
        O.elements[h].set_kf(k_hinges)
    for p in range(len(hinges), len(hinges)+len(panels)):
        O.elements[p].set_kf(k_panels)
    O.solve()
    O.exportJSON('./Examples/Mesh_tests/Miura_Unit_Cell.json')

    displacements = []
    load_factors = []
    for i in range(len(O.solver.solutions)):
        O.solver.setSolution(i, elements=True)
        displacements.append([-O.U[5*3][0]])
        load_factors.append([O.solution_info['ld']])
        print(
            f"Displacement: {-O.U[5*3][0]} Load factor: {O.solution_info['ld']}")
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

    html = anim.to_jshtml()
    with open('./mori.html', 'w') as f:

        f.write(html)
