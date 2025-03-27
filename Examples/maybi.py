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
    stifness_f = 500
    k_hinges = 1
    k_panels = stifness_f*k_hinges
    E = 1e4
    A = 1
    filename = "waterbomb.fold"
    input_data = BarAndHingeNonLinear.import_fold_file(filename)

    miura_coords = np.array(input_data[0])
    bars = np.array(input_data[1])
    hinges = np.array(input_data[2])
    mva = np.array(input_data[3])
    supports = [[0, 1, 1, 1],
                [172, 1, 1, 1],
                [15, 1, 1, 1]]

    hinge_bars = hinges[:, [1, 2]]

    elements = hinges.tolist() + hinge_bars.tolist()+bars.tolist()

    types = ['OH']*len(hinges) + ['L1V']*(len(hinge_bars)+len(bars))
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
    # for load in loads:
    # node, sx, sy, sz = load
    # O.addLoadNode(int(node), [float(sx), float(sy), float(sz)])
    for i, hinge in enumerate(hinges):
        try:
            force = -1
            if mva[i] == 0:
                force = 1
            O.elements[i].set_internal_force(force)
        except:
            pass

    O.solver.set_increments(100)
    O.solver.maxiter = 500
    O.solver.tol = 1e-3
    O.solver.set_delta_lambda_bar(10)
    fig = plt.figure()
    ax = fig.add_subplot()
    for i, e in enumerate(O.elements):
        if types[i] == "OH":
            # anotate hinge number in the middle of the hinge
            if i < len(hinges):
                x = np.mean(e.coords[:, 0])
                y = np.mean(e.coords[:, 1])
                # z = np.mean(e.coords[:, 2])
                if mva[i] == 0:
                    ax.text(x, y, f"{i}", color='red', fontsize=10)
                    ax.plot(e.coords[:, 0], e.coords[:, 1],
                            '--', c='r', zorder=10)
                else:
                    ax.text(x, y, f"{i}", color='blue', fontsize=10)
                    ax.plot(e.coords[:, 0], e.coords[:, 1],
                            '--', c='b', zorder=10)
        else:
            if i < len(hinges) or (i > len(hinges) and i < len(hinges)):
                ax.plot(e.coords[:, 0], e.coords[:, 1], 'k-')
    for i, node in enumerate(O.geometry.gdls):
        ax.text(node[0], node[1], f"{i}", fontsize=10)

    plt.axis('off')
    plt.show()

    for h in range(len(hinges)):
        O.elements[h].set_kf(k_hinges)
    for p in range(len(hinges), len(hinges)):
        O.elements[p].set_kf(k_panels)
    # O.solve()
    # O.exportJSON('./Examples/Mesh_tests/Fold_test_full_guaterbomb.json')

    solutions = json.load(
        open('./Examples/Mesh_tests/Fold_test_full_guaterbomb_300.json'))
    for i, sol in enumerate(solutions["solutions"]):
        U = np.array(sol["U"])
        info = sol["info"]
        if i != 16:
            O.solver.solutions.append(U)
            O.solver.solutions_info.append(info)
    O.solver.setSolution(-1, True)
    O.solve(guess=O.solver.solutions[-1], _guess=True)
    O.exportJSON('./Examples/Mesh_tests/Fold_test_full_guaterbomb.json')
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
        ax.set_box_aspect([1, 1, 1])
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
    ax.set(xlim3d=(0, 20), xlabel='X')
    ax.set(ylim3d=(0, 20), ylabel='Y')
    ax.set(zlim3d=(0, 20), zlabel='Z')
    # html = anim.to_jshtml()
    # with open('./MVassigment_unfolded_masgrande_miura2.html', 'w') as f:

    #     f.write(html)
