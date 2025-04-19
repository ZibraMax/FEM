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
    for fman in range(530, 1000, 10):

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

        O.solver.set_increments(10)
        O.solver.maxiter = 25
        O.solver.tol = 1e-3
        O.solver.set_delta_lambda_bar(0.125)

        for h in range(len(hinges)):
            O.elements[h].set_kf(k_hinges)
        for p in range(len(hinges), len(hinges)):
            O.elements[p].set_kf(k_panels)

        solutions = json.load(
            open(f'./Examples/Mesh_tests/Fold_test_full_guaterbomb_{fman}.json'))
        for i, sol in enumerate(solutions["solutions"]):
            U = np.array(sol["U"])
            info = sol["info"]
            if i != 16:
                O.solver.solutions.append(U)
                O.solver.solutions_info.append(info)
        O.solver.setSolution(-1, True)
        O.solve(guess=O.solver.solutions[-1], _guess=True)
        O.exportJSON(
            f'./Examples/Mesh_tests/Fold_test_full_guaterbomb_{fman+10}.json')
