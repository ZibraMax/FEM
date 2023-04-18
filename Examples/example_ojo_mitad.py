if __name__ == '__main__':
    from FEM.Geometry import Delaunay, Region2D
    from FEM.Elasticity2D import PlaneStrainSparse
    import numpy as np
    import matplotlib.pyplot as plt
    W = 3.5*2
    pata = 8
    r = 12
    th = 1
    r_int = r - th
    n = 50
    A = 0.5
    pressure = 100
    E = 1000
    v = 0.3
    ancho_nervio = 2.5

    essential1 = 49
    essential2 = 106

    def circulo(x):
        centro = pata+r
        x = x - centro
        if (r**2-x**2) < 0:
            return W/2
        return max((r**2-x**2)**0.5, W/2)

    def circulo_int(x):
        centro = pata+r
        x = x - centro
        if (r_int**2-x**2) < 0:
            return 0
        return (r_int**2-x**2)**0.5

    def contour_ext(x):
        return circulo(x)

    x_des = -(r_int**2-(ancho_nervio)**2)**0.5 + pata + r

    def contour_int(x):
        return circulo_int(x)

    def different_region(x):
        return ancho_nervio

    X = np.linspace(0, pata+r, n)
    _X = []
    _Y = []
    for x in X:
        _X.append(x)
        _Y.append(contour_ext(x))
    flag = True
    for x in X[::-1]:
        if x >= pata+th:
            if x < x_des and flag:
                flag = False
                _X.append(x_des)
                _Y.append(contour_int(x_des))
            else:
                _X.append(x)
                _Y.append(contour_int(x))
    if flag:
        _X.append(x_des)
        _Y.append(contour_int(x_des))
    _X.append(pata+th)
    _Y.append(0)

    vertices = []
    for x, y in zip(_X, _Y):
        vertices.append([x, y])
    vertices.append([0, -ancho_nervio])

    params = Delaunay._strdelaunay(True, True, A)
    o = Delaunay(vertices, params, 2, fast=True)
    o.show()
    plt.show()
    ngeo = o.revolve(m=30)
    ngeo.exportJSON("EYE_REVOLVED.json")
