"""Polygonal help functions
"""


import numpy as np
import random
import math
from typing import Tuple


def enmalladoFernando(lx: float, ly: float, nex: int, ney: int, filename: str) -> None:
    """Crea un enmallado 2D de un rectangulo

    Args:
            lx (floar): Base del rectángulo
            ly (float): Altura del rectámgulo
            nex (int): Numero de elementos en el eje x
            ney (int): Numero de elementos en el eje y
            filename (str): nombre del archivo donde se guardará
    """

    lx = float(lx)
    ly = float(ly)
    nex = int(nex)
    ney = int(ney)
    lr = 6*l
    nne = 8
    hx = lx/nex
    hy = ly/ney
    nnd = (ney+1)*(2*nex+1)+(ney)*(nex+1)
    x = np.zeros([nnd])
    y = np.zeros([nnd])
    nel = nex*ney
    nle = []
    elm = np.zeros([nel, 9])
    ntel = 1
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
        elm[ne, 8] = 1
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
            elm[ne, 8] = 1

    print('Guardando Archivo')
    f = open(filename, 'w')
    f.write(format(nnd)+'\t'+format(nel)+'\t0\t0\t0\t2'+'\n')
    for i in range(nnd):
        f.write(format(x[i])+'\t'+format(y[i])+'\n')
    for i in range(nel):
        f.write('C2V'+'\n')
    for i in range(nel):
        def fun(x): return str(int(x)-1)
        f.write('\t'.join(map(fun, [elm[i, 0], elm[i, 1], elm[i, 2],
                elm[i, 3], elm[i, 4], elm[i, 5], elm[i, 6], elm[i, 7]]))+'\n')
    f.close()
    print('Archivo ' + filename + ' Guardado')


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
    sum = 0
    for i in range(numVerts):
        tmp = random.uniform(lower, upper)
        angleSteps.append(tmp)
        sum = sum + tmp

    # normalize the steps so that point 0 and point n+1 are the same
    k = sum / (2*math.pi)
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


def clip(x: float, min: float, max: float) -> float:
    """Clip 1D

    Args:
        x (float): Point x
        min (float): min
        max (float): max

    Returns:
        float: idk
    """

    if(min > max):
        return x
    elif(x < min):
        return min
    elif(x > max):
        return max
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

    return math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)


def isBetween(a: list, b: list, c: list, tol: float) -> bool:
    """Test if a point is between a line in a given tolerance

    Args:
        a (list): Start point of line
        b (list): End point of line
        c (list): Point to be tested between line
        tol (float): Tolerance

    Returns:
        bool: True if point is in line
    """

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
    angle = (np.arctan2(P[1]-P1[1], P[0]-P1[0]) -
             np.arctan2(P[1]-P2[1], P[0]-P2[0]))/2
    segment = r/np.abs(np.tan(angle))
    PC1 = segment
    PC2 = segment
    PP1 = np.sqrt((P[0]-P1[0])**2+(P[1]-P1[1])**2)
    PP2 = np.sqrt((P[0]-P2[0])**2+(P[1]-P2[1])**2)
    minr = min(PP1, PP2)/2
    if segment > minr:
        segment = minr
        r = segment*np.abs(np.tan(angle))
    PO = np.sqrt(r**2+segment**2)
    C1 = [-1, -1]
    C2 = [-1, -1]
    C1[0] = P[0]-(P[0]-P1[0])*PC1/PP1
    C1[1] = P[1]-(P[1]-P1[1])*PC1/PP1

    C2[0] = P[0]-(P[0]-P2[0])*PC2/PP2
    C2[1] = P[1]-(P[1]-P2[1])*PC2/PP2

    C = [-1, -1]
    C[0] = C1[0]+C2[0]-P[0]
    C[1] = C1[1]+C2[1]-P[1]
    dx = P[0]-C[0]
    dy = P[1]-C[1]
    PC = np.sqrt(dx**2+dy**2)

    O = [-1, -1]
    O[0] = P[0]-dx*PO/PC
    O[1] = P[1]-dy*PO/PC

    sangle = np.arctan((C1[1]-O[1])/(C1[0]-O[0]))
    eangle = np.arctan((C2[1]-O[1])/(C2[0]-O[0]))

    sweep = eangle - sangle
    if sweep < 0.0:
        sweep = -sweep
        sangle = eangle
    if sweep > np.pi:
        sweep = np.pi-sweep
    return O, sangle, sweep


def giveCoordsCircle(O: list, r: float, sa: float = 0, a: float = np.pi*2, n: int = 10) -> Tuple[np.ndarray, list]:
    """Calculates the coordinates of a circle

    Args:
        O (list): Center coordinates of circle
        r (float): Circle radius
        sa (float): Start angle. Defaults to 0
        a (float): End angle. Defaults to :math: `2\\pi`
        n (int, optional): Number of coords to calculate. Defaults to 10.

    Returns:
        np.ndarray and list: Circle coordinates and segments
    """
    coords = []
    segments = []
    h = a/n
    for i in range(n):
        if i < n-1:
            segments += [[i, i+1]]
        else:
            segments += [[i, 0]]
        theta = sa+h*i
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        coords += [[O[0]+x, O[1]+y]]
    return np.array(coords), segments
