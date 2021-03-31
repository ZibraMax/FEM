"""General geometry class.
"""


import numpy as np
from ..Utils import isBetween, roundCorner, giveCoordsCircle, angleBetweenAngles
import matplotlib.pyplot as plt
from ..Elements.E1D.LinealElement import LinealElement
from ..Elements.E1D.QuadraticElement import QuadraticElement
from ..Elements.E2D.Serendipity import Serendipity
from ..Elements.E2D.Quadrilateral import Quadrilateral
from ..Elements.E2D.QTriangular import QTriangular
from ..Elements.E2D.LTriangular import LTriangular
from typing import Callable
from ast import literal_eval
from tqdm import tqdm
import re


class Geometry:
    """Define geometry structure

    Args:
        dictionary (list): Matrix with element definitions. Each row is an element. The gdl are defined in columns
        gdls (list): List of domain coordinates
        types (list): Types of each element
        nvn (int, optional): Nunmber of variables per node. Defaults to 1.
        segments (list, optional): Domain segments. Defaults to [].
    """

    def __init__(self, dictionary: list, gdls: list, types: list, nvn: int = 1, segments: list = []) -> None:
        """Define geometry structure

        Args:
            dictionary (list): Matrix with element definitions. Each row is an element. The gdl are defined in columns
            gdls (list): List of domain coordinates
            types (list): Types of each element
            nvn (int, optional): Nunmber of variables per node. Defaults to 1.
            segments (list, optional): Domain segments. Defaults to [].
        """

        self.mask = None
        self.holes = None
        self.fillets = None
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

    def maskFromSegments(self) -> None:
        """Create the display mask from geometry segments
        """

        self.mask = []
        for s in self.segments:
            self.mask += np.array(self.gdls)[np.ix_(s)].tolist()

    def initialize(self) -> None:
        """Calculates the total number of GDL's and generates the elements structure
        """

        self.ngdl = int(len(self.gdls)*self.nvn)
        self.generateElements()

    def detectNonLocal(self, lr: float) -> list:
        """Detect adjacent elements between a distance Lr

        Args:
            lr (float): Distance to detect adjacent elements

        Returns:
            list: Non local element dictionary
        """

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
    def loadGiDMsh(filename: str):
        """Load geometry from GiD msh file

        Args:
            filename (str): Path to GiD file

        Returns:
            Geometry: Result geometry
        """

        f = open(filename, 'r')
        dicc = []
        gdls = []
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
    def loadmsh(filename: str):
        """Load geometry from previously generated MSH file

        Args:
            filename (str): Path to msh file

        Returns:
            Geometry: Output geometry
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
        try:
            len_holes = p[6]
        except:
            len_holes = 0
        try:
            len_fillets = p[7]
        except:
            len_fillets = 0
        try:
            len_mask = p[8]
        except:
            len_mask = 0
        fillets = []
        holes = []
        mask = []
        for _ in range(len_holes):
            holes.append(literal_eval(f.readline()))
        for _ in range(len_fillets):
            fillets.append(literal_eval(f.readline()))
        for _ in range(len_mask):
            mask += [list(map(float, f.readline().split('\t')))]
        f.close()
        print('File ' + filename + ' loaded')
        o = Geometry(dicc, gdls, types, nvn, seg)
        o.cbe = cbe
        o.cbn = cbn
        o.holes = holes
        o.fillets = fillets
        o.mask = mask
        return o

    def generateElements(self) -> None:
        """Generate elements structure
        """
        print('Generating element structure')
        self.elements = []
        for i, d in enumerate(tqdm(self.dictionary, unit='Element')):
            coords = np.array(self.gdls)[np.ix_(d)]
            gdl = np.zeros([self.nvn, len(d)])
            for j in range(self.nvn):
                gdl[j, :] = (np.array(d)*self.nvn+j)
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
        print('Done!')

    def show(self, texto: int = 10, bolita: int = 0, draw_segs: bool = True, draw_labels: bool = False, draw_bc: bool = False, label_bc: bool = False) -> None:
        """Create a geometry graph

        Args:
            texto (int, optional): Text size. Defaults to 10.
            bolita (int, optional): Node size. Defaults to 0.
            draw_segs (bool, optional): To draw or not draw the segments. Defaults to True.
            draw_labels (bool, optional): To draw or not draw element labels. Defaults to False.
            draw_bc (bool, optional): To draw border conditions. Defaults to False.
            label_bc (bool, optional): To draw labels on border conditions. Defaults to False.
        """

        fig = plt.figure()
        ax = fig.add_subplot()

        ax.axes.set_aspect('equal')

        for i, e in enumerate(self.elements):
            coords = e._coords
            coords = np.array(coords.tolist() + [coords[0].tolist()])
            X = coords[:, 0]
            Y = coords[:, 1]
            ax.plot(X, Y, '-', color='black', alpha=1-0.6*draw_bc, zorder=-10)
            cx = self.centroids[i][0]
            cy = self.centroids[i][1]
            if draw_labels:
                ax.plot(cx, cy, 'o', markersize=texto + bolita, color='yellow')
                ax.annotate(format(i), [
                            cx, cy], size=texto, textcoords="offset points", xytext=(-0, -2.5), ha='center')
        try:
            if draw_segs:
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
                        alpha=1-0.6*draw_bc
                    )
                    cx = (x0+x1)*0.5
                    cy = (y0+y1)*0.5
                    ax.plot(cx, cy, 'o', markersize=texto +
                            bolita, color='pink', alpha=1-0.6*draw_bc)
                    ax.annotate(format(i), [
                                cx, cy], alpha=1-0.6*draw_bc, size=texto, textcoords="offset points", xytext=(-0, -2.5), ha='center')
        except:
            pass
        for j, e in enumerate(self.elements):
            if e.intBorders:
                for i in range(-1, len(e.borders)-1):
                    border = e.borders[i]
                    if (len(border.properties['load_x']) + len(border.properties['load_y'])):
                        coords_border_0 = e._coords[i]
                        coords_border_1 = e._coords[i+1]
                        ax.plot([coords_border_0[0], coords_border_1[0]], [coords_border_0[1], coords_border_1[1]],
                                color='yellow',
                                linewidth=5,
                                zorder=50,
                                )
                        cx = (coords_border_0 + coords_border_1)/2
                        ax.annotate(format(len(border.properties['load_x']) + len(border.properties['load_y'])), cx, size=texto,
                                    textcoords="offset points", xytext=(-0, -2.5), ha='center', zorder=55)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('Domain')

        gdls = np.array(self.gdls)

        labels = np.linspace(0, gdls.shape[0] - 1, gdls.shape[0]).astype(int)
        if draw_labels:
            ax.plot(gdls[:, 0], gdls[:, 1], 'o',
                    markersize=texto+bolita, color='gray')
        if draw_labels:
            for p, l in zip(gdls, labels):
                ax.annotate(l, p, size=texto, textcoords="offset points",
                            xytext=(-0, -2.5), ha='center')
        maxx = np.max(gdls[:, 0])
        maxy = np.max(gdls[:, 1])

        minx = np.min(gdls[:, 0])
        miny = np.min(gdls[:, 1])
        coordmax = min(maxx-minx, maxy-miny)
        tFlecha = coordmax/80
        if draw_bc:
            for i, cb in enumerate(self.cbe):
                coords_cb = gdls[int(cb[0]//self.nvn)]
                if cb[0] % self.nvn == 0:
                    color = 'red'
                    ax.annotate(f"{i}: {cb[1]}"*label_bc, xy=coords_cb, xytext=(
                        coords_cb[0]-tFlecha, coords_cb[1]), horizontalalignment='center', verticalalignment='center', arrowprops=dict(arrowstyle="->", facecolor=color))
                elif cb[0] % self.nvn == 1:
                    color = 'blue'
                    ax.annotate(f"{i}: {cb[1]}"*label_bc, xy=coords_cb, xytext=(
                        coords_cb[0], coords_cb[1]+tFlecha), horizontalalignment='center', verticalalignment='center', arrowprops=dict(arrowstyle="->", facecolor=color))
                elif cb[0] % self.nvn == 2:
                    color = 'yellow'
                    ax.annotate(f"{i}: {cb[1]}"*label_bc, xy=coords_cb, xytext=(
                        coords_cb[0]-tFlecha, coords_cb[1]-tFlecha), horizontalalignment='center', verticalalignment='center', arrowprops=dict(arrowstyle="->", facecolor=color))
                else:
                    color = 'black'
                    ax.annotate(f"{i}: {cb[1]}"*label_bc, xy=coords_cb, xytext=(
                        coords_cb[0]-tFlecha, coords_cb[1]), horizontalalignment='center', verticalalignment='center', arrowprops=dict(arrowstyle="->", facecolor=color))
        figManager = plt.get_current_fig_manager()
        figManager.full_screen_toggle()

    def saveMesh(self, ProjectName: str) -> None:
        """Saves the geometry to a MSH file with specified name

        Args:
            ProjectName (str): Project name without extension
        """
        filename = ProjectName + '.msh'
        f = open(filename, 'w')
        try:
            len_holes = len(self.holes)
        except Exception as e:
            len_holes = 0
        try:
            len_fillets = len(self.fillets)
        except Exception as e:
            len_fillets = 0
        try:
            len_mask = len(self.mask)
        except Exception as e:
            len_mask = 0
        p = [len(self.gdls), len(self.dictionary), len(
            self.segments), len(self.cbe), len(self.cbn), self.nvn, len_holes, len_fillets, len_mask]
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
        if self.holes:
            for e in self.holes:
                f.write(str(e)+'\n')
        if self.fillets:
            for e in self.fillets:
                f.write(str(e)+'\n')
        if self.mask:
            for e in self.mask:
                f.write('\t'.join(list(map(format, e)))+'\n')
        f.close()
        print('File ' + filename + ' saved')

    def centroidsAndAreas(self) -> None:
        """Calculate elements centroids and areas
        """
        for e in self.elements:
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

    def setCbe(self, cbe: list) -> None:
        """This method have to be used to assign essential boundary conditions. Thes method prevents to assign duplicated border conditions

        Args:
            cbe (list): Border conditions to be applied
        """
        res = []
        for i in cbe:
            if i not in res:
                res.append(i)
        self.cbe = res

    def generateSegmentsFromCoords(self, p0: list, p1: list) -> None:
        """Generates a geometry segment by specified coordinates

        Args:
            p0 (list): Segment start point
            p1 (list): Segment end point
        """
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

    def generateBCFromCoords(self, x: float, y: float, value: float = 0, nv: int = 1) -> list:
        """Generates border conditions by coordinates. The border condition is applied to the nearest node

        Args:
            x (float): X coordinate of point
            y (float): Y coordinate of point
            value (float, optional): Value of the border condition. Defaults to 0.
            nv (int, optional): Variable number. The first variable is 1. Defaults to 1.

        Returns:
            list: Matrix of border coordinates that can be concatenated
        """

        masCercano1 = None
        d1 = np.Inf
        for i, gdl in enumerate(self.gdls):
            r1 = np.sqrt((x-gdl[0])**2+(y-gdl[1])**2)
            if r1 < d1:
                d1 = r1
                masCercano1 = i
        return [[masCercano1*self.nvn+(nv-1), value]]

    def giveNodesOfSegment(self, segment: int, tol: float) -> np.ndarray:
        """Give nodes over a segment

        Args:
            segment (int): Segment number. Start with 0
            tol (float): Tolerance for finding nearest nodes

        Returns:
            np.ndarray: List of nodes in the specified segment
        """
        a = []
        ps = np.array(self.gdls)[self.segments[segment]].tolist()
        for i, p in enumerate(self.gdls):
            if isBetween(ps[0], ps[1], p, tol):
                a.append(i)
        return np.array(a)

    def giveElementsOfSegment(self, segment: int, tol: float) -> list:
        """Give elements over a segment

        Args:
            segment (int): Segment number. Start with 0
            tol (float): Tolerance for finding nearest nodes

        Returns:
            lsit: List of elements in the specified segment
        """
        a = []
        nodes = self.giveNodesOfSegment(segment, tol)
        for e in self.elements:
            if np.sum(np.isin(e.gdl[0], nodes*self.nvn)) > 0:
                a.append(e)
        return a

    def loadOnSegment(self, segment: int, fx: Callable = None, fy: Callable = None, tol: float = 1*10**(-5), add=None) -> None:
        """Assign a load over a geometry segment.

        The start point of segment is the 0 point of load
        The end point of segment is the end point of load

        Load must be defined as a function (normal or lambda)

        Args:
            segment (int): Segment in wich load will be applied
            fx (Callable, optional): Load Function x component. Defaults to None.
            fy (Callable, optional): Load Function y component. Defaults to None.
            tol (float, optional): Tolerancy for finding nodes. Defaults to 1*10**(-5).
        """
        a = self.giveElementsOfSegment(segment, tol)
        coordenadas = np.array(self.gdls)[self.segments[segment]]
        vect_seg = coordenadas[1]-coordenadas[0]
        for e in a:
            e.intBorders = True
            for i in range(-1, len(e.borders)-1):
                pertenece1 = isBetween(
                    coordenadas[0], coordenadas[1], e._coords[i], tol)
                pertenece2 = isBetween(
                    coordenadas[0], coordenadas[1], e._coords[i+1], tol)
                if pertenece1 and pertenece2:
                    vect_lad = e._coords[i+1]-e._coords[i]
                    sign = np.sign(vect_seg@vect_lad)
                    e.borders[i].dir = sign
                    e.borders[i].s0 = np.linalg.norm(
                        e._coords[i]-coordenadas[0])
                    if fx:
                        e.borders[i].properties['load_x'].append(fx)
                    if fy:
                        e.borders[i].properties['load_y'].append(fy)
                    if add:
                        e.borders[i].properties.update(add)
                else:
                    e.borders[i].dir = 0.0

    def loadOnHole(self, hole: int, sa: float = 0, ea: float = 2*np.pi, fx: Callable = None, fy: Callable = None, tol: float = 1*10**(-5)) -> None:
        """Assign loads over a hole.

        Args:
            hole (int): Hole index in wich load will be applied
            sa (float, optional): Start face angle. Defaults to 0.
            ea (float, optional): Finish face angle. Defaults to :math:`2\\pi`.
            fx (Callable, optional): Load Function x component. Defaults to None.
            fy (Callable, optional): Load Function y component. Defaults to None.
            tol (float, optional): Tolerancy for finding nodes. Defaults to 1*10**(-5).
        """
        holee = self.holes[hole]
        segments_apply = []
        for i, segmento in enumerate(holee['segments']):
            seg_coords = np.array(self.gdls)[segmento]
            centradas = seg_coords[1]-seg_coords[0]
            angle = np.arctan2(centradas[1], centradas[0])
            angle += np.pi/2
            if angle < 0:
                angle += 2*np.pi
            if angleBetweenAngles(sa, ea, angle):
                segments_apply.append(segmento)
        for segmento in segments_apply:
            for i, seg in enumerate(self.segments):
                if seg == segmento:
                    self.loadOnSegment(i, fx, fy, tol)
                    break

    def cbOnHole(self, hole: int, value: float, nv: int = 1, sa: float = 0, ea: float = 2*np.pi, tol: float = 1*10**(-5)) -> list:
        """Generate a list of border conditions from specified hole.

        Args:
            hole (int): Hole index in wich load will be applied
            value (float): Value of the bc
            nv (int, optional): Variable number, starts with 1. Defaults to 1.
            sa (float, optional): Start face angle. Defaults to 0.
            ea (float, optional): Finish face angle. Defaults to :math:`2\\pi`.
            tol (float, optional): Tolerancy for finding nodes. Defaults to 1*10**(-5).

        Returns:
            list: List of border conditions that can be concatenated or assigned to the geometry
        """
        holee = self.holes[hole]
        segments_apply = []
        bc = []
        for i, segmento in enumerate(holee['segments']):
            seg_coords = np.array(self.gdls)[segmento]
            centradas = seg_coords[1]-seg_coords[0]
            angle = np.arctan2(centradas[1], centradas[0])
            angle += np.pi/2
            if angle < 0:
                angle += 2*np.pi
            if angleBetweenAngles(sa, ea, angle):
                segments_apply.append(segmento)
        for segmento in segments_apply:
            for i, seg in enumerate(self.segments):
                if seg == segmento:
                    bc += self.cbFromSegment(i, value, nv, tol)
                    break
        return bc

    def cbFromSegment(self, segment: int, value: float, nv: int = 1, tol: float = 1*10**(-5)) -> list:
        """Generate a list of border conditions from specified border.

        Args:
            segment (int): Segment number
            value (float): Value of the bc
            nv (int, optional): Variable number, starts with 1. Defaults to 1.
            tol (float, optional): Tolerance for finding nodes in the segment. Defaults to 1*10**(-5).

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

    def cbeAllBorders(self, value: float, tol: float = 1*10**(-5)):
        """Set all segments border conditions to the specified value

        Args:
            value (float): Value of the border condition
            tol (float, optional): Tolerance for finding near nodes in segments. Defaults to 1*10**(-5).
        """
        for s in range(len(self.segments)):
            for i in range(self.nvn):
                self.cbe += self.cbFromSegment(s, value, (i+1), tol)
