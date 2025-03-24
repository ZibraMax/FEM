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
    E = 1e4
    A = 1

    input_data = json.load(open("input_data.json"))

    miura_coords = np.array(input_data["nodes"])
    perimeter = np.array(input_data["perimeter"]).astype(int)-1
    hinges = np.array(input_data["folds"]).astype(int)-1
    panels = np.array(input_data["bend"]).astype(int)-1

    supports = np.array(input_data["supports"]).astype(int)
    supports[:, 0] -= 1

    loads = np.array(input_data["loads"])
    loads[:, 0] -= 1

    hinges = hinges[:, [2, 0, 1, 3]]
    panels = panels[:, [2, 0, 1, 3]]

    hinge_bars = hinges[:, [1, 2]]
    panel_bars = panels[:, [1, 2]]

    elements = hinges.tolist() + panels.tolist() + perimeter.tolist() + \
        hinge_bars.tolist() + panel_bars.tolist()
    types = ['OH']*len(hinges) + ['OH']*len(panels) + ['L1V'] * \
        len(perimeter) + ['L1V']*len(hinge_bars) + ['L1V']*len(panel_bars)
    geo = Geometry3D(elements, miura_coords.tolist(), types, 3, fast=True)
    for e in supports:
        node, sx, sy, sz = e
        if sx == 1:
            geo.cbe.append([int(node)*3, 0.0])
        if sy == 1:
            geo.cbe.append([int(node)*3+1, 0.0])
        if sz == 1:
            geo.cbe.append([int(node)*3+2, 0.0])
    O = BarAndHingeNonLinear(geo, E, A, verbose=True)
    for load in loads:
        node, sx, sy, sz = load
        O.addLoadNode(int(node), [float(sx), float(sy), float(sz)])
    O.solver.set_increments(60)
    O.solver.maxiter = 100
    O.solver.tol = 1e-6
    O.solver.set_delta_lambda_bar(0.5)
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    for i, e in enumerate(O.elements):
        if types[i] == "OH":
            ax.plot(e.coords[:, 0], e.coords[:, 1],
                    e.coords[:, 2], '--', c='yellow', zorder=10)
        else:
            ax.plot(e.coords[:, 0], e.coords[:, 1],
                    e.coords[:, 2], 'k-')
    plt.show()

    for h in range(len(hinges)):
        O.elements[h].set_kf(k_hinges)
    for p in range(len(hinges), len(hinges)+len(panels)):
        O.elements[p].set_kf(k_panels)
    O.solve()
    O.exportJSON('./Examples/Mesh_tests/Bar_and_hinge_please_work.json')

    displacements = []
    load_factors = []
    for i in range(len(O.solver.solutions)):
        O.solver.setSolution(i, elements=True)
        O.elements[0].calculate_vectors()
        theta = O.elements[0].calculate_theta()
        displacements.append([theta])
        load_factors.append([O.solution_info['ld']])
        print(
            f"Angle: {theta} Load factor: {O.solution_info['ld']}")
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
    html = anim.to_jshtml()
    with open('./M_URO.html', 'w') as f:

        f.write(html)
