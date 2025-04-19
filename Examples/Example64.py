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

    import sympy as sp

    beta, a, b, gamma1, gamma2 = sp.symbols(
        "\\beta a b \\gamma_1 \\gamma_2", real=True, positive=True)

    psi = sp.acos(-sp.cos(gamma1)*sp.cos(gamma2) +
                  sp.sin(gamma1)*sp.sin(gamma2)*sp.cos(beta))
    exp = sp.cos(beta)-(sp.sin(beta)**2*sp.sin(gamma1)
                        * sp.sin(gamma2))/(1-sp.cos(psi))
    theta = sp.acos(exp)
    theta = sp.simplify(theta)
    phi = sp.acos(sp.cos(gamma1)*sp.cos(gamma2) +
                  sp.sin(gamma1)*sp.sin(gamma2)*sp.cos(theta))
    phi = sp.simplify(phi)
    L = a*sp.sin(psi/2)
    m = a*sp.cos(psi/2)
    W = b*sp.sin(phi/2)
    V = b*sp.cos(phi/2)
    c1 = a*sp.sin(gamma2)
    c2 = a*sp.sin(gamma1)
    LL = sp.sqrt(c1**2+c2**2-2*c1*c2*sp.cos(beta))
    k1 = sp.asin(sp.sin(beta)*c2/LL)
    H = c1*sp.sin(k1)

    LL = sp.sqrt(c1**2+c2**2-2*c1*c2*sp.cos(beta))

    fH = sp.lambdify([beta, gamma1, gamma2, a, b], H, cse=True)
    fL = sp.lambdify([beta, gamma1, gamma2, a, b], 2*L, cse=True)
    fV = sp.lambdify([beta, gamma1, gamma2, a, b], V, cse=True)
    fS = sp.lambdify([beta, gamma1, gamma2, a, b], 2*W, cse=True)
    ftheta = sp.lambdify([beta, gamma1, gamma2, a, b], theta, cse=True)

    funcs = [fH, fL, fV, fS]

    def R(theta):
        return theta*np.pi/180

    def D(theta):
        return theta*180/np.pi
    stifness_f = 999999999
    k_hinges = 1
    k_panels = stifness_f*k_hinges
    E = 999999999
    A = 1

    input_data = json.load(open("input_data_unfolded.json"))

    miura_coords = np.array(input_data["nodes"])
    perimeter = np.array(input_data["perimeter"]).astype(int)-1
    hinges = np.array(input_data["folds"]).astype(int)-1
    panels = np.array(input_data["bend"]).astype(int)-1

    supports = np.array(input_data["supports"]).astype(int)
    supports[:, 0] -= 1
    assigment = input_data["assigment"]

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
    # for load in loads:
    # node, sx, sy, sz = load
    # O.addLoadNode(int(node), [float(sx), float(sy), float(sz)])
    for i, hinge in enumerate(hinges):
        force = -1
        if assigment[i] == "M":
            force = 1
        O.elements[i].set_internal_force(force)

    O.solver.set_increments(60)
    O.solver.maxiter = 100
    O.solver.tol = 1e-6
    O.solver.set_delta_lambda_bar(0.5)
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    for i, e in enumerate(O.geometry.gdls):
        ax.plot(e[0], e[1], e[2], 'o', markersize=5, c='k')
        ax.text(e[0], e[1], e[2],
                f"{i}",
                size=10, zorder=11, color='red', ha='center', va='center')

    for i, e in enumerate(O.elements):
        if types[i] == "OH":
            ax.plot(e.coords[:, 0], e.coords[:, 1],
                    e.coords[:, 2], '--', c='yellow', zorder=10)
        else:
            ax.plot(e.coords[:, 0], e.coords[:, 1],
                    e.coords[:, 2], 'k-')
            center = e._xcenter
            # annotate in the center
            ax.text(center[0], center[1], center[2],
                    f"{i}",
                    size=10, zorder=11, color='blue', ha='center', va='center')
    plt.show()

    for h in range(len(hinges)):
        O.elements[h].set_kf(k_hinges)
    for p in range(len(hinges), len(hinges)+len(panels)):
        O.elements[p].set_kf(k_panels)
    O.solve()
    O.exportJSON('./Examples/Mesh_tests/Fold_test_2.json')
    properties = {}
    properties["H"] = []
    properties["V"] = []
    properties["L"] = []
    properties["W"] = []
    displacements = []
    load_factors = []
    for i in range(len(O.solver.solutions)):
        O.solver.setSolution(i, elements=True)
        lit = []
        for k in range(4):
            O.elements[k].calculate_vectors()
            theta = O.elements[k].calculate_theta()*180/np.pi
            lit.append(theta)
        V = O.geometry.gdls[3][1] + O.U[10][0]
        H = O.U[5][0]
        L = O.geometry.gdls[2][1] + O.U[7][0]
        W = O.geometry.gdls[6][0] + O.U[18][0]
        properties["V"].append(V)
        properties["H"].append(H)
        properties["L"].append(L)
        properties["W"].append(W)
        displacements.append(lit)
        load_factors.append([O.solution_info['ld']]*len(lit))
        print(
            f"Angle: {displacements[-1]} Load factor: {O.solution_info['ld']}")
    displacements = np.array(displacements)
    load_factors = np.array(load_factors)
    # Plot the results
    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig.add_subplot(1, 2, 2)
    cosas = ax2.plot(displacements, load_factors, 'o-',
                     label=[f"$\\rho_{i}$" for i in range(len(displacements[0]))])

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
        ax.set_box_aspect([1, 1, 1])
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
    with open('./MVassigment_2.html', 'w') as f:
        f.write(html)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    i = 0
    ts = np.linspace(0, np.pi, 100)
    for prop, values in properties.items():
        ax.plot(displacements[:, 1], values, '-', label=prop, lw=3, zorder=10)
        thetas_para_miura_eq = ftheta(ts, R(60), R(60), 2, 2)
        aa = funcs[i](ts, R(60), R(60), 2, 2)
        ax.plot(thetas_para_miura_eq*180/np.pi, aa,
                '--', label=f"{prop} (analitical)", lw=3, zorder=30)
        i += 1
    ax.set_xlabel('Angle')
    ax.set_ylabel('Displacement')
    ax.legend()
    ax.grid()
    plt.show()
