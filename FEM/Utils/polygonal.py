import math
import random
import numpy as np


def enmalladoFernando(l, lx, ly, nex, ney, filename):
    """Crea un enmallado 2D de un rectangulo

    Args:
            l (float): Parámetro no local no usado en este programa, dejar en 0
            lx (floar): Base del rectángulo
            ly (float): Altura del rectámgulo
            nex (int): Numero de elementos en el eje x
            ney (int): Numero de elementos en el eje y
            filename (str): nombre del archivo donde se guardará
    """
    l = float(l)
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
    # Identify neighbors within influence distance
    print('Detectando Elementos No Locales')
    nem = 0
    for i in range(0, nel):
        ne = -1
        nii = elm[i, 0]
        nfi = elm[i, 2]
        xci = (x[int(nii-1)]+x[int(nfi-1)])/2
        yci = (y[int(nii-1)]+y[int(nfi-1)])/2
        ne = ne+1
        g = ne
        nnn = []
        for j in range(0, nel):
            if not j == i:
                nij = elm[j, 0]
                nfj = elm[j, 2]
                xcj = (x[int(nij-1)]+x[int(nfj-1)])/2
                ycj = (y[int(nij-1)]+y[int(nfj-1)])/2
                dist = np.sqrt((xci-xcj)**2+(yci-ycj)**2)
                if dist < (lr+0.5*np.sqrt(hx**2+hy**2)):
                    ne = ne+1
                    nnn.append(j+1)
        nle.append(np.array([ne]+[i+1]+nnn))
        if ne > nem:
            nem = ne
    nle = np.array(nle)
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


def generatePolygon(ctrX=10, ctrY=10, aveRadius=5, irregularity=0.5, spikeyness=0.5, numVerts=6):
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


def clip(x, min, max):
    if(min > max):
        return x
    elif(x < min):
        return min
    elif(x > max):
        return max
    else:
        return x


def dist(a, b):
    return math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)


def isBetween(a, b, c, tol):
    d1 = dist(a, c)
    d2 = dist(b, c)
    d3 = dist(a, b)
    d = d1+d2-d3
    if abs(d) < tol:
        return True
    return False
