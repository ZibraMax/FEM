if __name__ == '__main__':
    from FEM.Geometry import Delaunay, Region2D
    from FEM.Elasticity2D import PlaneStrainSparse
    import numpy as np
    import matplotlib.pyplot as plt
    W = 4
    pata = 3
    r = 12
    th = 1
    r_int = r - th
    n = 50
    A = 0.1
    pressure = 100
    E = 21000000
    v = 0.3

    essential1 = 49
    essential2 = 124

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

    x_des = -(r_int**2-(W/4)**2)**0.5 + pata + r

    def contour_int(x):
        return circulo_int(x)

    def different_region(x):
        return W/4

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
    for x, y in zip(_X[::-1][1:], _Y[::-1][1:]):
        vertices.append([x, -y])

    vertices.append([0, -W/4])
    vertices.append([0, W/4])
    n = len(vertices)
    for i, p in enumerate(vertices):
        if p[0] == x_des and p[1] == contour_int(x_des):
            e1 = i
        elif p[0] == x_des and p[1] == -contour_int(x_des):
            e2 = i

    s1 = [n-1, e1]
    s2 = [n-2, e2]

    es = [s1, s2]

    reg = Region2D(np.array([
        [0.0, W/4],
        [pata+th, W/4],
        [pata+th, -W/4],
        [0.0, -W/4]]))

    params = Delaunay._strdelaunay(True, True, A, o=2)
    geo = Delaunay(vertices, params, 2, extra_segs=es, fast=True)
    geo.addRegions([reg])
    cb = geo.cbFromRegion(essential1, 0.0, 1)
    cb += geo.cbFromRegion(essential2, 0.0, 1)
    cb += geo.generateBCFromCoords(pata+th, 0, 0, 2)
    geo.setCbe(cb)
    for i in range(essential1+1, essential2):
        geo.normalLoadOnRegionVF(i, lambda s: -pressure)
    # geo.show()
    # plt.show()
    for e in geo.giveElementsOfRegion(-1, True):
        e.properties['E'] = 0.1*E
        e.properties['v'] = v
    ES = []
    VS = []
    for e in geo.elements:
        ES.append(e.properties.get("E", E))
        VS.append(e.properties.get("v", v))

    o = PlaneStrainSparse(geo, ES, VS, verbose=True)
    o.solve()
    o.exportJSON("EYE_1.json")
    plt.show()
