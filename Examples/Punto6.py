if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from FEM.Geometry.Geometry import Geometry
    from FEM.Poisson2D import Poisson2D
    phi = 2
    nodos = [[0.0, 0.0], [0.5, 0.0], [1.0, 0.0],
             [0.0, 0.5], [0.5, 0.5], [1.0, 0.5],
             [0.5, 1.0], [1.0, 1.0], [1.0, 1.5]]
    elementos = [[0, 1, 4, 3], [1, 2, 5, 4],
                 [3, 4, 6], [4, 5, 7, 6], [6, 7, 8]]
    tipos = ['C1V', 'C1V', 'T1V', 'C1V', 'T1V']
    geometry = Geometry(elementos, nodos, tipos)
    geometry.cbe = [[2, 0], [5, 0], [7, 0], [8, 0], [0, 18], [3, 18]]
    geometry.cbn = [[6, 0], [1, 0]]
    # geometry.show(draw_segs=True, draw_labels=True, draw_bc=True, label_bc=True)
    # plt.show()
    O = Poisson2D(geometry, phi)
    O.solve(plot=False)
    # plt.show()
    plt.show()
