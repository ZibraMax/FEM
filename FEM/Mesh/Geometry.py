"""This should be documented
"""


import numpy as np
from ..Utils import isBetween
import matplotlib.pyplot as plt
from ..Elements import *
from ..Elements.E1D import *
from ..Elements.E2D import *
# from ..Elements.E3D import *
import re


class Geometry:
    """Aver si esto funciona aqui."""

    def __init__(self, dictionary, gdls, types, nvn=1, segments=[]):
        """Aver si esto funciona aqui.

        Args:
            dictionary ([type]): [description]
            gdls ([type]): [description]
            types ([type]): [description]
            nvn (int, optional): [description]. Defaults to 1.
            segments (list, optional): [description]. Defaults to [].
        """
        self.mask = None
        self.areas = []
        self.nvn = nvn
        self.dictionary = dictionary
        self.elements = []
        self.gdls = gdls
        self.types = types
        self.segments = segments
        self.cbe = []
        self.cbn = []
        self.centroids = []
        self.initialize()
        try:
            self.centroidsAndAreas()
        except:
            pass

    def maskFromSegments(self):
        self.mask = []
        for s in self.segments:
            self.mask += np.array(self.gdls)[np.ix_(s)].tolist()

    def initialize(self):
        self.ngdl = int(len(self.gdls)*self.nvn)
        self.generateElements()

    def detectNonLocal(self, lr):
        diccionariosnl = []
        for i in range(len(self.dictionary)):
            cxl, cyl = self.centroids[i]
            linea = []
            linea.append(i)
            for j in range(len(self.dictionary)):
                if not j == i:
                    cxnl, cynl = self.centroids[j]
                    d = ((cxl-cxnl)**2+(cyl-cynl)**2)**0.5
                    if d <= lr:
                        linea.append(j)
            diccionariosnl.append(linea)
        return diccionariosnl

    @staticmethod
    def loadGiDMsh(filename):
        f = open(filename, 'r')
        dicc = []
        gdls = []
        types = []
        seg = []
        cbe = []
        cbn = []
        nvn = 1
        params = re.findall(r'\S+', f.readline())
        f.readline()
        while True:
            string = f.readline()
            if string.split('\n')[0] == 'End Coordinates':
                break
            gdls += [list(map(float, re.findall(r'\S+', string)[1:3]))]
        f.readline()
        f.readline()
        while True:
            string = f.readline()
            if string.split('\n')[0] == 'End Elements':
                break
            dicc += [list(map(lambda x: int(x)-1,
                          re.findall(r'\S+', string)[1:]))]
        f.close()
        tipo = 'C2V'
        if params[4] == 'Quadrilateral':
            tipo = 'C1V'
            if params[-1] == '8':
                tipo = 'C2V'
        if params[4] == 'Triangle':
            tipo = 'T1V'
            if params[-1] == '6':
                tipo = 'T2V'
        types = [tipo]*len(dicc)
        print('File ' + filename + ' loaded')
        o = Geometry(dicc, gdls, types)
        return o

    @staticmethod
    def loadmsh(filename):
        """Aver, esto debería funcionar re pro

        Args:
            filename ([type]): [description]

        Returns:
            [type]: [description]
        """
        print('Loading ' + filename)
        f = open(filename, 'r')
        dicc = []
        gdls = []
        types = []
        seg = []
        cbe = []
        cbn = []
        nvn = 1
        p = list(map(int, re.findall(r'\S+', f.readline())))
        # [len(self.gdls),len(self.dictionary),len(self.segments),len(self.cbe),len(self.cbn),self.nvn]
        for _ in range(p[0]):
            gdls += [list(map(float, f.readline().split('\t')))]
        for _ in range(p[1]):
            types += [f.readline().split('\n')[0]]
        for _ in range(p[1]):
            dicc += [list(map(int, f.readline().split('\t')))]
        for _ in range(p[2]):
            seg += [list(map(int, f.readline().split('\t')))]
        for _ in range(p[3]):
            cbe += [list(map(float, f.readline().split('\t')))]
        for _ in range(p[4]):
            cbn += [list(map(float, f.readline().split('\t')))]
        nvn = p[5]
        f.close()
        print('File ' + filename + ' loaded')
        o = Geometry(dicc, gdls, types, nvn, seg)
        o.cbe = cbe
        o.cbn = cbn
        return o

    def generateElements(self):
        self.elements = []
        for i, d in enumerate(self.dictionary):
            coords = np.array(self.gdls)[np.ix_(d)]
            gdl = np.zeros([self.nvn, len(d)])
            for i in range(self.nvn):
                gdl[i, :] = (np.array(d)*self.nvn+i)
            gdl = gdl.astype(int)
            if self.types[i] == 'T1V':
                element = LTriangular(coords, gdl)
            elif self.types[i] == 'T2V':
                element = QTriangular(coords, gdl)
            elif self.types[i] == 'C1V':
                element = Quadrilateral(coords, gdl)
            elif self.types[i] == 'C2V':
                element = Serendipity(coords, gdl)
            elif self.types[i] == 'L1V':
                element = LinealElement(coords, gdl)
            elif self.types[i] == 'L2V':
                element = QuadraticElement(coords, gdl)
            self.elements.append(element)

    def show(self, texto=10, bolita=0, figsize=[17, 10]):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot()

        ax.axes.set_aspect('equal')

        for i, e in enumerate(self.elements):
            coords = e._coords
            coords = np.array(coords.tolist() + [coords[0].tolist()])
            X = coords[:, 0]
            Y = coords[:, 1]
            ax.plot(X, Y, 'o-', color='black', zorder=-10)
            cx = self.centroids[i][0]
            cy = self.centroids[i][1]
            ax.plot(cx, cy, 'o', markersize=texto + bolita, color='yellow')
            ax.annotate(format(i), [
                        cx, cy], size=texto, textcoords="offset points", xytext=(-0, -2.5), ha='center')
        try:
            verts = self.gdls
            segs = self.segments
            for i, seg in enumerate(segs):
                x0, y0 = verts[int(seg[0])]
                x1, y1 = verts[int(seg[1])]

                ax.fill(
                    [x0, x1],
                    [y0, y1],
                    facecolor='none',
                    edgecolor='b',
                    linewidth=3,
                    zorder=0,
                )
                cx = (x0+x1)*0.5
                cy = (y0+y1)*0.5
                ax.plot(cx, cy, 'o', markersize=texto + bolita, color='pink')
                ax.annotate(format(i), [
                            cx, cy], size=texto, textcoords="offset points", xytext=(-0, -2.5), ha='center')
        except:
            pass
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('Domain')

        gdls = np.array(self.gdls)

        labels = np.linspace(0, gdls.shape[0] - 1, gdls.shape[0]).astype(int)

        ax.plot(gdls[:, 0], gdls[:, 1], 'o',
                markersize=texto+bolita, color='gray')

        for p, l in zip(gdls, labels):
            ax.annotate(l, p, size=texto, textcoords="offset points",
                        xytext=(-0, -2.5), ha='center')

    def saveMesh(self, ProjectName):
        filename = ProjectName + '.msh'
        f = open(filename, 'w')
        p = [len(self.gdls), len(self.dictionary), len(
            self.segments), len(self.cbe), len(self.cbn), self.nvn]
        f.write('\t'.join(list(map(format, p))) + '\n')
        for e in self.gdls:
            f.write('\t'.join(list(map(format, e))) + '\n')
        for e in self.types:
            f.write(e + '\n')
        for e in self.dictionary:
            f.write('\t'.join(list(map(format, e))) + '\n')
        for e in self.segments:
            f.write('\t'.join(list(map(format, e))) + '\n')
        for e in self.cbe:
            f.write(format(int(e[0]))+'\t'+format(e[1])+'\n')
        for e in self.cbn:
            f.write(format(int(e[0]))+'\t'+format(e[1])+'\n')
        f.close()
        print('File ' + filename + ' saved')

    def centroidsAndAreas(self):
        for i, e in enumerate(self.elements):
            coords = e._coords
            coords = np.array(coords.tolist() + [coords[0].tolist()])
            area = 0
            cx = 0
            cy = 0
            for j in range(len(coords)-1):
                area += coords[j][0]*coords[j+1][1]-coords[j+1][0]*coords[j][1]
                mult = (coords[j][0]*coords[j+1][1] -
                        coords[j+1][0]*coords[j][1])
                cx += (coords[j][0]+coords[j+1][0])*mult
                cy += (coords[j][1]+coords[j+1][1])*mult
            self.areas.append(np.abs(area/2))
            self.centroids.append([cx/3/area, cy/3/area])

    def generateSegmentsFromCoords(self, p0, p1):
        masCercano1 = None
        d1 = np.Inf
        masCercano2 = None
        d2 = np.Inf
        for i, gdl in enumerate(self.gdls):
            r1 = np.sqrt((p0[0]-gdl[0])**2+(p0[1]-gdl[1])**2)
            r2 = np.sqrt((p1[0]-gdl[0])**2+(p1[1]-gdl[1])**2)
            if r1 < d1:
                d1 = r1
                masCercano1 = i
            if r2 < d2:
                d2 = r2
                masCercano2 = i
        self.segments.append([masCercano1, masCercano2])

    def generateBCFromCoords(self, x, y, value=0, nv=1):
        masCercano1 = None
        d1 = np.Inf
        for i, gdl in enumerate(self.gdls):
            r1 = np.sqrt((x-gdl[0])**2+(y-gdl[1])**2)
            if r1 < d1:
                d1 = r1
                masCercano1 = i
        return [[masCercano1*self.nvn+(nv-1), value]]

    def giveNodesOfSegment(self, segment, tol):
        a = []
        ps = np.array(self.gdls)[self.segments[segment]].tolist()
        for i, p in enumerate(self.gdls):
            if isBetween(ps[0], ps[1], p, tol):
                a.append(i)
        return np.array(a)

    def cbFromSegment(self, segment, value, nv=1, tol=1*10**(-5)):
        """Generate a list of border conditions from specified border.

        Args:
            segment (int): Segment number
            value (float): Value of the bc
            nv (int, optional): Variable number, starts with 1. Defaults to 1.
            tol (float, optional): Tolerancy for finding nodes in the segment. Defaults to 1*10**(-5).

        Returns:
            list: List of border conditions that can be concatenated or assigned to the geometry
        """
        cb = []
        nodes = self.giveNodesOfSegment(segment, tol)
        cbe = np.zeros([len(nodes), 2])
        cbe[:, 0] = nodes*self.nvn+(nv-1)
        cbe[:, 1] = value
        cb += cbe.tolist()
        return cb

    def cbeAllBorders(self, value, tol=1*10**(-5)):
        for s in range(len(self.segments)):
            for i in range(self.nvn):
                self.cbe += self.cbFromSegment(s, value, (i+1), tol)
