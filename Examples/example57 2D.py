if __name__ == '__main__':
    from FEM import TrussNonLinear, TrussLinear
    from FEM.Geometry import Geometry3D
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.animation import FuncAnimation

    a = 1
    b = 1
    L = (a**2 + b**2)**0.5
    h = 0.5
    t = 0.1
    A = h*t
    P = 1000
    young = 210e3  # Young's modulus in MPa
    EA = young*A

    coords = [
        [0, 0, 0.0],
        [a, b, 0.0],
        [2*a, 0, 0.0]
    ]
    elements = [[0, 1],
                [1, 2]]
    types = ['L1V']*(len(elements))

    geo = Geometry3D(elements, coords, types, 3, fast=True)
    for node in [0, 2]:
        geo.cbe += [[node*3, 0], [node*3+1, 0], [node*3+2, 0]]
    geo.cbe += [[1*3+2, 0]]
    O = TrussNonLinear(geo, EA, 1)
    O.solver.set_delta_lambda_bar(0.1)
    O.solver.set_increments(220)
    O.solver.momentum = False
    O.addLoadNode(1, [0.0, -1000, 0])
    O.solve()
    O.exportJSON('./Examples/Mesh_tests/Truss_non_lineal2D.json')

    displacements = []
    load_factors = []
    for i in range(len(O.solver.solutions)):
        O.solver.setSolution(i, elements=True)
        displacements.append(-O.U[4][0])
        load_factors.append(O.solution_info['ld'])

    data = np.array([displacements, load_factors]).T
    np.savetxt('./Examples/examples_results/bar_and_hinge_formulation.csv', data,
               header='Displacement\tLoad factor', delimiter=',')

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
