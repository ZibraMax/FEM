if __name__ == '__main__':
    from FEM.Heat2D import Heat2D
    from FEM.Geometry import Delaunay
    import matplotlib.pyplot as plt
    # https://github.com/Samson-Mano/2D_Heat_transfer
    beta1 = 2
    beta2 = 0.8
    k = 0.8

    ta1 = 10
    ta2 = 150

    b = 80
    h = 100

    b_in = 50
    h_in = 70

    tx = (b-b_in)/2
    ty = (h-h_in)/2

    coords = [[0, 0], [b, 0], [b, h], [0, h]]
    coords_hueco = [[tx, ty], [b-tx, ty], [b-tx, h-ty], [tx, h-ty]]
    seg_hueco = [[0, 1], [1, 2], [2, 3], [3, 0]]
    params = Delaunay._strdelaunay(constrained=True, delaunay=True,
                                   a='5', o=2)
    cent = [b/2, h/2]
    hole = {'center': cent, 'regions': seg_hueco, 'vertices': coords_hueco}
    holes = [hole]

    geometria = Delaunay(coords, params, 1, holes_dict=holes)
    O = Heat2D(geometria, k, k)
    O.defineConvectiveBoderConditions(0, beta2, ta1)
    O.defineConvectiveBoderConditions(1, beta2, ta1)
    O.defineConvectiveBoderConditions(2, beta2, ta1)
    O.defineConvectiveBoderConditions(3, beta2, ta1)

    O.defineConvectiveBoderConditions(4, beta1, ta2)
    O.defineConvectiveBoderConditions(5, beta1, ta2)
    O.defineConvectiveBoderConditions(6, beta1, ta2)
    O.defineConvectiveBoderConditions(7, beta1, ta2)
    O.solve()
    O.exportJSON("Examples/Mesh_tests/Example36.json")
    plt.show()
