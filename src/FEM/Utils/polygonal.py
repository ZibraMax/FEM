"""Polygonal help functions
"""


import numpy as np
import matplotlib.pyplot as plt
import random
import math
from typing import Tuple


def enmalladoEsferaFernando(L: float, n: float) -> Tuple[np.ndarray, list]:
    """Crea el enmallado de una esfera de diámetro L con n numero de elementos
        por lado. Para crear el enmallado se deforma un cubo a la forma de una
        esfera.

        Autor: Fernando Ramirez Rodriguez
        Traducido de Matlab.

    Args:
        L (float): Diametro de la esfera
        n (float): Número de elementos por lado. El numero deelemtnos final es n^3

    Returns:
        Tuple[np.ndarray,list]: Matriz de coordenadas y conectividad
    """
    # TODO Necesita revisión

    sidel = L
    nsdel = n
    delta = sidel/nsdel

    re = sidel/2
    deltar = 2*re/nsdel
    nd = 0
    coor = np.zeros([(nsdel+1)**3, 3])
    for i in range(0, nsdel+1):
        z = (i)*delta
        for j in range(0, nsdel+1):
            y = (j)*delta
            for k in range(0, nsdel+1):
                x = (k)*delta
                coor[nd, 0] = x-sidel/2
                coor[nd, 1] = y-sidel/2
                coor[nd, 2] = z-sidel/2
                nd = nd+1
    ntnd = nd
    tol = 1e-6
    for i in range(0, int(nsdel/2)):
        crl = sidel/2-(i)*delta
        rec = re-(i)*deltar
        for nd in range(0, ntnd):
            if abs(abs(coor[nd, 0]) - crl) < tol or abs(abs(coor[nd, 1]) - crl) < tol or abs(abs(coor[nd, 2]) - crl) < tol:
                d = np.sqrt((coor[nd, 0]) ** 2+(coor[nd, 1])
                            ** 2+(coor[nd, 2]) ** 2)
                dd = rec-d
                xu = coor[nd, 0]/d
                yu = coor[nd, 1]/d
                zu = coor[nd, 2]/d
                coor[nd, 0] = coor[nd, 0]+xu*dd
                coor[nd, 1] = coor[nd, 1]+yu*dd
                coor[nd, 2] = coor[nd, 2]+zu*dd
    for nd in range(0, ntnd):
        coor[nd, 1-1] = coor[nd, 1-1]+sidel/2
        coor[nd, 2-1] = coor[nd, 2-1]+sidel/2
        coor[nd, 3-1] = coor[nd, 3-1]+sidel/2
    # return coor
    con = []
    el = 0
    for i in range(nsdel):
        for j in range(nsdel):
            for k in range(nsdel):
                ni = (i)*(nsdel+1)*(nsdel+1)+(j)*(nsdel+1)+k
                con.append([0]*8)
                con[el][1-1] = ni
                con[el][2-1] = ni+1
                con[el][3-1] = con[el][2-1]+nsdel+1
                con[el][4-1] = con[el][1-1]+nsdel+1
                con[el][5-1] = con[el][1-1]+(nsdel+1)*(nsdel+1)
                con[el][6-1] = con[el][5-1]+1
                con[el][7-1] = con[el][6-1]+nsdel+1
                con[el][8-1] = con[el][5-1]+nsdel+1
                el = el + 1
    return coor, con


def enmalladoFernando(lx: float, ly: float, nex: int, ney: int) -> np.ndarray:
    """Crea un enmallado 2D de un rectangulo

    Args:
            lx (floar): Base del rectángulo
            ly (float): Altura del rectámgulo
            nex (int): Numero de elementos en el eje x
            ney (int): Numero de elementos en el eje y
    Returns:
            np.ndarray: coordinates matrix (np.ndarray) and element dictionary (list)
    """

    lx = float(lx)
    ly = float(ly)
    nex = int(nex)
    ney = int(ney)
    hx = lx/nex
    hy = ly/ney
    nnd = (ney+1)*(2*nex+1)+(ney)*(nex+1)
    x = np.zeros([nnd])
    y = np.zeros([nnd])
    nel = nex*ney
    elm = np.zeros([nel, 8])
    # Coordinate Generation
    print('Generando Coordenadas')
    nd = -1
    for i in range(1, ney+1):
        cy = (i-1)*hy
        for j in range(1, 2*nex+2):
            nd = nd+1
            y[nd] = cy
            x[nd] = (j-1)*hx/2
        cy = (i-1)*hy+hy/2
        for j in range(1, nex+2):
            nd = nd+1
            y[nd] = cy
            x[nd] = (j-1)*hx
    cy = ly
    for j in range(1, 2*nex+2):
        nd = nd+1
        y[nd] = cy
        x[nd] = (j-1)*hx/2
    # Element Node Connectivty
    print('Generando Elementos')

    ne = -1
    for i in range(0, ney):
        ne = ne+1
        elm[ne, 0] = (i)*(3*nex+2)+1
        elm[ne, 1] = elm[ne, 0]+2
        elm[ne, 3] = elm[ne, 0]+3*nex+2
        elm[ne, 2] = elm[ne, 3]+2
        elm[ne, 4] = elm[ne, 0]+1
        elm[ne, 7] = elm[ne, 0]+2*nex+1
        elm[ne, 6] = elm[ne, 3]+1
        elm[ne, 5] = elm[ne, 7]+1
        for j in range(1, nex):
            ne = ne+1
            elm[ne, 0] = elm[ne-1, 1]
            elm[ne, 1] = elm[ne, 0]+2
            elm[ne, 3] = elm[ne-1, 2]
            elm[ne, 2] = elm[ne, 3]+2
            elm[ne, 4] = elm[ne, 0]+1
            elm[ne, 7] = elm[ne-1, 5]
            elm[ne, 6] = elm[ne, 3]+1
            elm[ne, 5] = elm[ne, 7]+1

    # print('Guardando Archivo')
    # f = open(filename, 'w')
    # f.write(format(nnd)+'\t'+format(nel)+'\t0\t0\t0\t2'+'\n')
    # for i in range(nnd):
    #     f.write(format(x[i])+'\t'+format(y[i])+'\n')
    # for i in range(nel):
    #     f.write('C2V'+'\n')
    # for i in range(nel):
    #     def fun(x): return str(int(x)-1)
    #     f.write('\t'.join(map(fun, [elm[i, 0], elm[i, 1], elm[i, 2],
    #             elm[i, 3], elm[i, 4], elm[i, 5], elm[i, 6], elm[i, 7]]))+'\n')
    # f.close()
    coords = np.array([x, y]).T
    dicc = (elm-1).astype(int).tolist()
    return coords, dicc


def generatePolygon(ctrX: float = 10, ctrY: float = 10, aveRadius: float = 5, irregularity: float = 0.5, spikeyness: float = 0.5, numVerts: float = 6) -> list:
    """Generate a random polygon.

    Args:
        ctrX (float, optional): X centroid. Defaults to 10.
        ctrY (float, optional): Y centroid. Defaults to 10.
        aveRadius (float, optional): Average radious. Defaults to 5.
        irregularity (float, optional): Irregularity. Defaults to 0.5.
        spikeyness (float, optional): Spikeyness. Defaults to 0.5.
        numVerts (float, optional): Number of vertices. Defaults to 6.

    Returns:
        list: Poligon coordinates matrix.
    """

    irregularity = clip(irregularity, 0, 1) * 2*math.pi / numVerts
    spikeyness = clip(spikeyness, 0, 1) * aveRadius

    # generate n angle steps
    angleSteps = []
    lower = (2*math.pi / numVerts) - irregularity
    upper = (2*math.pi / numVerts) + irregularity
    suma = 0
    for i in range(numVerts):
        tmp = random.uniform(lower, upper)
        angleSteps.append(tmp)
        suma = suma + tmp

    # normalize the steps so that point 0 and point n+1 are the same
    k = suma / (2*math.pi)
    for i in range(numVerts):
        angleSteps[i] = angleSteps[i] / k

    # now generate the points
    points = []
    angle = random.uniform(0, 2*math.pi)
    for i in range(numVerts):
        r_i = clip(random.gauss(aveRadius, spikeyness), 0, 2*aveRadius)
        x = ctrX + r_i*math.cos(angle)
        y = ctrY + r_i*math.sin(angle)
        points.append((int(x), int(y)))

        angle = angle + angleSteps[i]

    return points


def clip(x: float, mi: float, ma: float) -> float:
    """Clip 1D

    Args:
        x (float): Point x
        mi (float): min
        ma (float): max

    Returns:
        float: idk
    """

    if (mi > ma):
        return x
    elif (x < mi):
        return mi
    elif (x > ma):
        return ma
    else:
        return x


def dist(a: list, b: list) -> float:
    """Calculate the distancie between 2 points

    Args:
        a (list): point a
        b (list): point b

    Returns:
        float: Distance between a and b
    """
    return np.linalg.norm(np.array(a)-np.array(b))


def isBetween(a: list, b: list, c: list, tol: float = 1*10**-5) -> bool:
    """Test if a point is between a line in a given tolerance. Works in 2D and 3D.

    Args:
        a (list): Start point of line
        b (list): End point of line
        c (list): Point to be tested between line
        tol (float): Tolerance. Defaults to 1*10**-5

    Returns:
        bool: True if point is in line
    """
    a = a.flatten()
    b = b.flatten()
    c = c.flatten()

    d1 = dist(a, c)
    d2 = dist(b, c)
    d3 = dist(a, b)
    d = d1+d2-d3
    if abs(d) < tol:
        return True
    return False


def roundCorner(P1: list, P2: list, P: list, r: float) -> tuple:
    """Calculates the origin, start angle and sweep angle of a given corner with a given radius
    Source: https://stackoverflow.com/questions/24771828/algorithm-for-creating-rounded-corners-in-a-polygon

    Args:
        P1 (list): First point
        P2 (list): Second Point
        P (list): Center point
        r (float): Radius of corner

    Returns:
        tuple: Circle center coordinates, start angle, sweep angle
    """
    def GetProportionPoint(point, segment, length, dx, dy):
        factor = segment / length
        return [point[0] - dx * factor, point[1] - dy * factor]
    dx1 = P[0]-P1[0]
    dy1 = P[1]-P1[1]
    dx2 = P[0]-P2[0]
    dy2 = P[1]-P2[1]

    angle = (np.arctan2(dy1, dx1)-np.arctan2(dy2, dx2))/2
    tan = np.abs(np.tan(angle))
    segment = r/tan

    len1 = np.sqrt(dx1**2+dy1**2)
    len2 = np.sqrt(dx2**2+dy2**2)
    length = np.min([len1, len2])
    if segment > length:
        print('The fillet radius is big')
    p1Cross = GetProportionPoint(P, segment, len1, dx1, dy1)
    p2Cross = GetProportionPoint(P, segment, len2, dx2, dy2)

    dx = P[0]*2-p1Cross[0]-p2Cross[0]
    dy = P[1]*2-p1Cross[1]-p2Cross[1]

    L = (dx**2+dy**2)**0.5
    d = (segment**2+r**2)**0.5
    circlePoint = GetProportionPoint(P, d, L, dx, dy)
    sa = np.arctan2(p1Cross[1]-circlePoint[1], p1Cross[0]-circlePoint[0])
    ea = np.arctan2(p2Cross[1]-circlePoint[1], p2Cross[0]-circlePoint[0])
    s = ea-sa
    # if s < 0:
    #     sa = ea
    #     s = -s
    if s > np.pi:
        s = np.pi-s
    return circlePoint, sa, s


def giveCoordsCircle(O: list, r: float, sa: float = 0, a: float = np.pi*2, n: int = 10, isFillet: bool = False) -> Tuple[list, list]:
    """Calculates the coordinates of a circle

    Args:
        O (list): Center coordinates of circle
        r (float): Circle radius
        sa (float): Start angle. Defaults to 0
        a (float): End angle. Defaults to :math:`2\\pi`
        n (int, optional): Number of coords to calculate. Defaults to 10.
        isFillet (bool, optional): If the circle will be used as fillet. Defaults to False.

    Returns:
        list and list: Circle coordinates and regions
    """
    coords = []
    regions = []
    h = a/n
    if isFillet:
        for i in range(n+1):
            regions += [[i, i+1]]
            theta = sa+h*i
            x = r*np.cos(theta)
            y = r*np.sin(theta)
            coords += [[O[0]+x, O[1]+y]]
        theta = a
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        coords += [[O[0]+x, O[1]+y]]

    else:
        for i in range(n):
            if i < n-1:
                regions += [[i, i+1]]
            else:
                regions += [[i, 0]]
            theta = sa+h*i
            x = r*np.cos(theta)
            y = r*np.sin(theta)
            coords += [[O[0]+x, O[1]+y]]
    return coords, regions


def angleBetweenAngles(start: float, end: float, mid: float) -> bool:
    """Evaluates if a angle is between 2 angles

    Args:
        start (float): Start angle
        end (float): End angle
        mid (float): Angle to be evaluated

    Returns:
        bool: Tru if mid is between start-end
    """
    end = end - start + 2*np.pi if (end - start) < 0.0 else end - start
    mid = mid - start + 2*np.pi if (mid - start) < 0.0 else mid - start
    return (mid < end)


def testNeighborg(e1, e2):
    # Este es el número de vertices mínimos para que un elemento sea vecino de otro
    MIN_VERTICES = 3
    en_comun = 0
    for c in e2.coords:
        test = any(np.equal(e1.coords, c).all(1))
        if test:
            en_comun += 1
            if en_comun >= MIN_VERTICES:
                return True
    return False


def plot_list_elements(l, c="k", acum=False):
    if not acum:
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
    else:
        ax = plt.gca()
    for e in l:
        ax.plot(*e._xcenter, "o", c=c)
