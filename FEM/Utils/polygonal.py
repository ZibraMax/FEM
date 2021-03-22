"""Polygonal help functions
"""


import numpy as np
import random
import math


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
