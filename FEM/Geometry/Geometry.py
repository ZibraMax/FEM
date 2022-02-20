"""Geometry definitions.
"""

import triangle as tr
import copy
import numpy as np
import json
from ..Utils import isBetween, roundCorner, giveCoordsCircle, angleBetweenAngles
import matplotlib.pyplot as plt
from ..Elements.E1D.LinealElement import LinealElement
from ..Elements.E1D.CubicElement import CubicElement
from ..Elements.E1D.QuadraticElement import QuadraticElement
from ..Elements.E2D.Serendipity import Serendipity
from ..Elements.E2D.Quadrilateral import Quadrilateral
from ..Elements.E2D.QTriangular import QTriangular
from ..Elements.E2D.LTriangular import LTriangular
from ..Elements.E3D.Brick import Brick
from .Region import Region, Region1D, Region2D
from typing import Callable
from ast import literal_eval, parse
from tqdm import tqdm
import re


class Geometry:
    """Define a general geometry structure

    Args:
        dictionary (list): Matrix with element definitions. Each row is an element. The gdl are defined in columns
        gdls (list): List of domain coordinates
        types (list): Types of each element
        nvn (int, optional): Nunmber of variables per node. Defaults to 1.
        regions (list, optional): List of domain regions. Defaults to [].
        fast (bool): If True, the created elements will have have the fast propertie (see Element class docs)
    """

    def __init__(self, dictionary: list, gdls: list, types: list, nvn: int = 1, regions: list[Region] = None, fast=False) -> None:
        """Define geometry structure

        Args:
            dictionary (list): Matrix with element definitions. Each row is an element. The gdl are defined in columns
            gdls (list): List of domain coordinates
            types (list): Types of each element
            nvn (int, optional): Nunmber of variables per node. Defaults to 1.
            regions (list, optional): List of domain regions. Defaults to None.
            fast (bool): If True, the created elements will have have the fast propertie (see Element class docs)
        """

        self.mask = None
        self.holes = []
        self.fillets = []
        self.nvn = nvn
        self.dictionary = dictionary
        self.elements = []
        self.gdls = np.array(gdls)
        self.types = types
        self.regions = regions or []
        self.cbe = []
        self.cbn = []
        self.centroids = []
        self.fast = fast
        self.initialize()
        self.calculateCentroids()

    def calculateRegions(self) -> None:
        """Calculates the nodes of the geometry regions
        """

        for region in tqdm(self.regions, unit="Region"):
            region.setNodesOfRegion(self)

    def maskFromRegions(self) -> None:
        """Create the display mask from geometry regions
        """
        pass  # TODO this have to be moved to the Geometry2D class
        # self.mask = []
        # for s in self.regions:
        #     self.mask += np.array(self.gdls)[np.ix_(s)].tolist()

    def initialize(self) -> None:
        """Calculates the total number of GDL's and generates the elements structure
        """

        self.ngdl = int(len(self.gdls)*self.nvn)
        self.generateElements()
        self.calculateRegions()

    def detectNonLocal(self, lr: float) -> list:
        """Detect adjacent elements between a distance Lr

        Args:
            lr (float): Distance to detect adjacent elements

        Returns:
            list: Non local element dictionary
        """
        print('Detecting non local elements')
        # FIXME hacerlo mas eficiente.
        diccionariosnl = []
        for i in tqdm(range(len(self.dictionary)), unit='Elements'):
            ci = self.centroids[i]
            linea = []
            linea.append(i)
            for j in range(len(self.dictionary)):
                if not j == i:
                    cnl = self.centroids[j]
                    d = np.linalg.norm(np.array(cnl)-np.array(ci))
                    if d <= lr:
                        linea.append(j)
            diccionariosnl.append(linea)
        return diccionariosnl

    def generateElements(self) -> None:
        """Generate elements structure
        """
        print('Generating element structure')
        self.elements = [0.0]*len(self.dictionary)
        for i, d in enumerate(tqdm(self.dictionary, unit='Element')):
            coords = self.gdls[np.ix_(d)]
            gdl = np.zeros([self.nvn, len(d)])
            for j in range(self.nvn):
                gdl[j, :] = (np.array(d)*self.nvn+j)
            gdl = gdl.astype(int)
            if self.types[i] == 'T1V':
                element = LTriangular(coords, gdl, fast=self.fast)
            elif self.types[i] == 'T2V':
                element = QTriangular(coords, gdl, fast=self.fast)
            elif self.types[i] == 'C1V':
                element = Quadrilateral(coords, gdl, fast=self.fast)
            elif self.types[i] == 'C2V':
                element = Serendipity(coords, gdl, fast=self.fast)
            elif self.types[i] == 'L1V':
                element = LinealElement(coords, gdl, fast=self.fast)
            elif self.types[i] == 'L2V':
                element = QuadraticElement(coords, gdl, fast=self.fast)
            elif self.types[i] == 'L3V':
                element = CubicElement(coords, gdl, fast=self.fast)
            elif self.types[i] == 'B1V':
                element = Brick(coords, gdl, fast=self.fast)
            self.elements[i] = element
        print('Done!')

    def show(self) -> None:
        """Creates a geometry graph"""
        pass

    def calculateCentroids(self) -> None:
        """Calculate elements centroids
        """
        for e in self.elements:
            x, p = e.T(e.center.T)
            self.centroids.append(x.tolist())

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

    def giveNodesOfRegion(self, region: int) -> np.ndarray:
        """Give nodes over a region

        Args:
            region (int): region number. Start with 0

        Returns:
            np.ndarray: List of nodes in the specified region
        """
        return self.regions[region].nodes

    def giveElementsOfRegion(self, region: int) -> list:
        """Give elements over a region

        Args:
            region (int): region number. Start with 0

        Returns:
            list: List of elements in the specified region
        """
        a = []
        nodes = self.giveNodesOfRegion(region)
        for e in self.elements:
            if np.sum(np.isin(e.gdl[0], nodes*self.nvn)) > 0:
                a.append(e)
        return a

    def cbFromRegion(self, region: int, value: float, nv: int = 1) -> list:
        """Generate a list of border conditions from specified border.

        Args:
            region (int): region number
            value (float): Value of the bc
            nv (int, optional): Variable number, starts with 1. Defaults to 1.

        Returns:
            list: List of border conditions that can be concatenated or assigned to the geometry
        """

        cb = []
        nodes = self.giveNodesOfRegion(region)
        cbe = np.zeros([len(nodes), 2])
        cbe[:, 0] = nodes*self.nvn+(nv-1)
        cbe[:, 1] = value
        cb += cbe.tolist()
        return cb

    def cbeAllRegions(self, value: float) -> None:
        """Set all regions border conditions to the specified value to all the variables.

        Args:
            value (float): Value of the border condition
        """
        self.cbe = []
        for s in range(len(self.regions)):
            for i in range(self.nvn):
                self.cbe += self.cbFromRegion(s, value, (i+1))

    def exportJSON(self, filename: str = None) -> str:
        """Export geometry definition as JSON file or JSON string

        Args:
            filename (str, optional): If given, a JSON file is created. Defaults to None.

        Returns:
            str: JSON string
        """
        x = {
            "nodes": self.gdls.tolist(),
            "dictionary": self.dictionary,
            "types": self.types,
            "regions": self.giveRegions(),
            "ebc": self.cbe,
            "nbc": self.cbn,
            "nvn": self.nvn,
            "ngdl": self.ngdl,
            "holes": self.holes,
            "fillets": self.fillets,
        }
        y = json.dumps(x)
        if filename:
            with open(filename, "w") as f:
                f.write(y)
        return y

    def giveRegions(self) -> list:
        """Returns a list of regions coordinates matrix

        Returns:
            list: List of regions coordinates matrix
        """
        coords = []
        for reg in self.regions:
            coords += [reg.coords.tolist()]
        return coords

    def addRegions(self, regions: list[Region]) -> None:
        """Adds regions to an already created geometry

        Args:
            regions (list[Region]): Regions to be created
        """
        for r in tqdm(regions, unit='Regions'):
            r.setNodesOfRegion(self)
        self.regions += regions
        # self.calculateRegions()

    @classmethod
    def importJSON(self, filename: str, **kargs) -> 'Geometry':
        """Import geometry definition from JSON file

        Args:
            filename (str): Path to the JSON file

        Returns:
            Geometry: Geometry generated using the JSON file
        """
        with open(filename) as f:
            parsed = json.loads(f.read())
            dcc = parsed['dictionary']
            nodes = parsed['nodes']
            types = parsed['types']
            nvn = parsed['nvn']
            regions = []
            regions_parsed = parsed['regions']
            for coords in regions_parsed:
                if len(coords) == 2:
                    regions.append(Region1D(coords))
                elif len(coords) == 4:
                    regions.append(Region2D(coords))
            o = self(dcc, nodes, types, nvn, regions, **kargs)
            o.cbe = parsed['ebc']
            o.cbn = parsed['nbc']
            o.holes = parsed['holes']
            o.fillets = parsed['fillets']
            return o


class Geometry1D(Geometry):
    """Define an 1D geometry structure

    Args:
        dictionary (list): Matrix with element definitions. Each row is an element. The gdl are defined in columns
        gdls (list): List of domain coordinates
        types (list): Types of each element
        nvn (int, optional): Nunmber of variables per node. Defaults to 1.
        fast (bool, optional): If True, the created elements will have have the fast propertie (see Element class docs)
    """

    def __init__(self, dictionary: list, gdls: list, types: list, nvn: int = 1, fast=False) -> None:
        """Define an 1D geometry structure

        Args:
            dictionary (list): Matrix with element definitions. Each row is an element. The gdl are defined in columns
            gdls (list): List of domain coordinates
            types (list): Types of each element
            nvn (int, optional): Nunmber of variables per node. Defaults to 1.
            fast (bool, optional): If True, the created elements will have have the fast propertie (see Element class docs)
        """

        Geometry.__init__(self, dictionary, gdls, types, nvn, [], fast)

    def generateElements(self) -> None:
        """Generate elements structure
        """
        for i, d in enumerate(self.dictionary):
            coords = np.array(self.gdls)[np.ix_(d)]
            gdl = np.zeros([self.nvn, len(d)])
            for i in range(self.nvn):
                gdl[i, :] = (np.array(d)*self.nvn+i)
            gdl = gdl.astype(int)
            if self.types[i] == 'L1V':
                element = LinealElement(coords, gdl)
            elif self.types[i] == 'L2V':
                element = QuadraticElement(coords, gdl)
            self.elements.append(element)

    def show(self) -> None:
        """Create a geometry graph
        """
        pass


class Geometry2D(Geometry):
    """Creates a 2D geometry

        Args:
            dictionary (list): Matrix with element definitions. Each row is an element. The gdl are defined in columns
            gdls (list): List of domain coordinates
            types (list): Types of each element
            nvn (int, optional): Nunmber of variables per node. Defaults to 1.
            regions (list[Region], optional): List of regions to apply in the geometry. Defaults to None.
            fast (bool, optional): If True, the created elements will have have the fast propertie (see Element class docs)
        """

    def __init__(self, dictionary: list, gdls: list, types: list, nvn: int = 1, regions: list[Region] = None, fast=False) -> None:
        """Creates a 2D geometry

        Args:
            dictionary (list): Matrix with element definitions. Each row is an element. The gdl are defined in columns
            gdls (list): List of domain coordinates
            types (list): Types of each element
            nvn (int, optional): Nunmber of variables per node. Defaults to 1.
            regions (list[Region], optional): List of regions to apply in the geometry. Defaults to None.
            fast (bool, optional): If True, the created elements will have have the fast propertie (see Element class docs)
        """
        Geometry.__init__(self, dictionary, gdls, types, nvn, regions, fast)

    def generateRegionFromCoords(self, p0: list, p1: list) -> None:
        """Generates a geometry Region1D by specified coordinates

        Args:
            p0 (list): region start point
            p1 (list): region end point
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
        coords = np.array([self.gdls[masCercano1], self.gdls[masCercano2]])
        self.regions.append(Region1D(coords))
        self.regions[-1].setNodesOfRegion(self)

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

    def loadOnRegionVF(self, region: int, f: Callable = None, add=None) -> None:
        """Assign a load over a geometry region.

        The start point of region is the 0 point of load
        The end point of region is the end point of load

        Load must be defined as a function (normal or lambda)

        Args:
            region (int): region in wich load will be applied
            f (Callable, optional): Load Function. Defaults to None.
        """
        c0, cf = self.regions[region].coords
        dy = cf[1]-c0[1]
        dx = cf[0]-c0[0]
        theta = np.arctan2(dy, dx)
        def fx(s): return f(c0[0]+s*np.cos(theta))[0]
        def fy(s): return f(c0[1]+s*np.sin(theta))[1]
        self.loadOnRegion(region=region, fx=fx,
                          fy=fy, add=add)

    def loadOnRegion(self, region: int, fx: Callable = None, fy: Callable = None, add=None) -> None:
        """Assign a load over a geometry region.

        The start point of region is the 0 point of load
        The end point of region is the end point of load

        Load must be defined as a function (normal or lambda)

        Args:
            region (int): region in wich load will be applied
            fx (Callable, optional): Load Function x component. Defaults to None.
            fy (Callable, optional): Load Function y component. Defaults to None.
        """
        a = self.giveElementsOfRegion(region)
        coordenadas = self.regions[region].coords
        vect_seg = coordenadas[1]-coordenadas[0]
        for e in a:
            e.intBorders = True
            for i in range(-1, len(e.borders)-1):
                pertenece1 = isBetween(
                    coordenadas[0], coordenadas[1], e._coords[i])
                pertenece2 = isBetween(
                    coordenadas[0], coordenadas[1], e._coords[i+1])
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

    def loadOnHole(self, hole: int, sa: float = 0, ea: float = 2*np.pi, fx: Callable = None, fy: Callable = None) -> None:
        """Assign loads over a hole.

        Args:
            hole (int): Hole index in wich load will be applied
            sa (float, optional): Start face angle. Defaults to 0.
            ea (float, optional): Finish face angle. Defaults to :math:`2\\pi`.
            fx (Callable, optional): Load Function x component. Defaults to None.
            fy (Callable, optional): Load Function y component. Defaults to None.
        """
        holee = self.holes[hole]
        regions_apply = []
        for i, region in enumerate(holee['regions']):
            seg_coords = self.gdls[region]
            centradas = seg_coords[1]-seg_coords[0]
            angle = np.arctan2(centradas[1], centradas[0])
            angle += np.pi/2
            if angle < 0:
                angle += 2*np.pi
            if angleBetweenAngles(sa, ea, angle):
                regions_apply.append(region)
        for region in regions_apply:
            for i, seg in enumerate(self.regions):
                if (seg.coords == self.gdls[np.ix_(region)]).all():
                    self.loadOnRegion(i, fx, fy)
                    break

    def cbOnHole(self, hole: int, value: float, nv: int = 1, sa: float = 0, ea: float = 2*np.pi) -> list:
        """Generate a list of border conditions from specified hole.

        Args:
            hole (int): Hole index in wich load will be applied
            value (float): Value of the bc
            nv (int, optional): Variable number, starts with 1. Defaults to 1.
            sa (float, optional): Start face angle. Defaults to 0.
            ea (float, optional): Finish face angle. Defaults to :math:`2\\pi`.

        Returns:
            list: List of border conditions that can be concatenated or assigned to the geometry
        """
        holee = self.holes[hole]
        regions_apply = []
        bc = []
        for i, region in enumerate(holee['regions']):
            seg_coords = self.gdls[region]
            centradas = seg_coords[1]-seg_coords[0]
            angle = np.arctan2(centradas[1], centradas[0])
            angle += np.pi/2
            if angle < 0:
                angle += 2*np.pi
            if angleBetweenAngles(sa, ea, angle):
                regions_apply.append(region)
        for region in regions_apply:
            for i, seg in enumerate(self.regions):
                if (seg.coords == self.gdls[np.ix_(region)]).all():
                    bc += self.cbFromRegion(i, value, nv)
                    break
        return bc

    def show(self, texto: int = 10, bolita: int = 0, draw_segs: bool = True, draw_labels: bool = False, draw_bc: bool = False, label_bc: bool = False) -> None:
        """Create a geometry graph

        Args:
            texto (int, optional): Text size. Defaults to 10.
            bolita (int, optional): Node size. Defaults to 0.
            draw_segs (bool, optional): To draw or not draw the regions. Defaults to True.
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
            ax.plot(X, Y, '-', color='black',
                    alpha=1-0.6*draw_bc, zorder=-10)
            cx = self.centroids[i][0][0]
            cy = self.centroids[i][0][1]
            if draw_labels:
                ax.plot(cx, cy, 'o', markersize=texto +
                        bolita, color='yellow')
                ax.annotate(format(i), [
                            cx, cy], size=texto, textcoords="offset points", xytext=(-0, -2.5), ha='center')
        try:
            if draw_segs:
                verts = self.gdls
                segs = self.regions
                for i, seg in enumerate(segs):
                    x0, y0 = seg.coords[0]
                    x1, y1 = seg.coords[1]

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
            for i, cb in enumerate(self.cbn):
                coords_cb = gdls[int(cb[0]//self.nvn)]
                if cb[0] % self.nvn == 0:
                    color = 'red'
                    ax.annotate(f"NBC {i}: {cb[1]}"*label_bc, xy=coords_cb, xytext=(
                        coords_cb[0]-tFlecha, coords_cb[1]), horizontalalignment='center', verticalalignment='center', arrowprops=dict(arrowstyle="->", facecolor=color))
                elif cb[0] % self.nvn == 1:
                    color = 'blue'
                    ax.annotate(f"NBC {i}: {cb[1]}"*label_bc, xy=coords_cb, xytext=(
                        coords_cb[0], coords_cb[1]+tFlecha), horizontalalignment='center', verticalalignment='center', arrowprops=dict(arrowstyle="->", facecolor=color))
                elif cb[0] % self.nvn == 2:
                    color = 'yellow'
                    ax.annotate(f"NBC {i}: {cb[1]}"*label_bc, xy=coords_cb, xytext=(
                        coords_cb[0]-tFlecha, coords_cb[1]-tFlecha), horizontalalignment='center', verticalalignment='center', arrowprops=dict(arrowstyle="->", facecolor=color))
                else:
                    color = 'black'
                    ax.annotate(f"NBC {i}: {cb[1]}"*label_bc, xy=coords_cb, xytext=(
                        coords_cb[0]-tFlecha, coords_cb[1]), horizontalalignment='center', verticalalignment='center', arrowprops=dict(arrowstyle="->", facecolor=color))
        figManager = plt.get_current_fig_manager()
        figManager.full_screen_toggle()


class Geometry3D(Geometry):
    """Creates a 2D geometry

    Args:
        dictionary (list): Matrix with element definitions. Each row is an element. The gdl are defined in columns
        gdls (list): List of domain coordinates
        types (list): Types of each element
        nvn (int, optional): Nunmber of variables per node. Defaults to 1.
        regions (list[Region], optional): List of regions to apply in the geometry. Defaults to None.
        fast (bool, optional): If True, the created elements will have have the fast propertie (see Element class docs)
    """

    def __init__(self, dictionary: list, gdls: list, types: list, nvn: int = 1, regions: list[Region] = None, fast=False):
        """Creates a 2D geometry

        Args:
            dictionary (list): Matrix with element definitions. Each row is an element. The gdl are defined in columns
            gdls (list): List of domain coordinates
            types (list): Types of each element
            nvn (int, optional): Nunmber of variables per node. Defaults to 1.
            regions (list[Region], optional): List of regions to apply in the geometry. Defaults to None.
            fast (bool, optional): If True, the created elements will have have the fast propertie (see Element class docs)
        """
        Geometry.__init__(self, dictionary, gdls, types, nvn, regions, fast)

    def show(self) -> None:
        """Creates a geometry graph
        """
        pass


class Lineal(Geometry1D):

    """Generate a evenly spaced elements domain

    Args:
        lenght (float): Domain lenght
        n (int): Number of elements
        o (int): Element order, can be 1 or 2
        nvn (int, optional): Number of variables per node. Defaults to 1.
    """

    def __init__(self, lenght: float, n: int, o: int, nvn: int = 1) -> None:
        """Generate a evenly spaced elements domain

        Args:
            lenght (float): Domain lenght
            n (int): Number of elements
            o (int): Element order, can be 1 or 2
            nvn (int, optional): Number of variables per node. Defaults to 1.
        """

        self.lenght = lenght
        dictionary = []
        gdls = []
        he = self.lenght / (n)
        for i in range(0, n):
            xa = i * he
            if o == 1:
                gdls += [xa]
                dictionary += [[i, i+1]]
            elif o == 2:
                gdls += [xa, xa+he/2]
                dictionary += [[i*o, i*o+1, i*o+2]]
            else:
                gdls += [xa, xa+he/3, xa+2*he/3]
                dictionary += [[i*o, i*o+1, i*o+2, i*o+3]]
        gdls += [self.lenght]
        if o == 1:
            tipo = 'L1V'
        elif o == 2:
            tipo = 'L2V'
        else:
            tipo = 'L3V'
        types = [tipo]*len(dictionary)
        gdls = np.array(gdls).reshape([len(gdls), 1])
        Geometry1D.__init__(self, dictionary, gdls, types,
                            nvn=nvn)


class Delaunay(Geometry2D):

    """Generate Delaunay triangulation using Triangle

    Args:
            vertices (list): matrix containing the domain vertices coordinates
            params (str): Triangulation parameters, use the aux function _strdelaunay
            nvn (int, optional): Number of variables per node. Defaults to 1.
            holes_dict (list, optional): A list of holes dicts. Defaults to None.
            fillets (list, optional): A list of fillets. Defaults to None.
            fast (bool, optional): If True, the created elements will have have the fast propertie (see Element class docs)
    """

    def __init__(self, vertices: list, params: str, nvn: int = 1, holes_dict=None, fillets=None, fast=False) -> None:
        """Generate Delaunay triangulation

        Args:
                vertices (list): matrix containing the domain vertices coordinates
                params (str): Triangulation parameters, use the aux function _strdelaunay
                nvn (int, optional): Number of variables per node. Defaults to 1.
                holes_dict (list, optional): A list of holes dicts. Defaults to None.
                fillets (list, optional): A list of fillets. Defaults to None.
                fast (bool, optional): If True, the created elements will have have the fast propertie (see Element class docs)
        """
        mask = copy.deepcopy(vertices)
        # try:
        # mask = mask.tolist()
        # except:
        # pass
        seg = []
        for i in range(len(vertices)-1):
            seg.append([i, i+1])
        seg.append([i+1, 0])
        hh = []
        mascarita = copy.deepcopy(seg)
        if fillets:
            for fillet in fillets:
                S1 = seg[fillet['start_region']]
                S2 = seg[fillet['end_region']]
                for i, maskarita in enumerate(mascarita):
                    if maskarita == S1:
                        indice_importante = i

                mizq = mascarita[:indice_importante]
                mder = mascarita[indice_importante:]

                P1 = vertices[S1[0]]
                P2 = vertices[S2[1]]
                P = vertices[S1[1]]
                r = fillet['r']
                n = fillet['n']
                if not S1[1] == S2[0]:
                    raise Exception('The fillet regions are not valid')
                O, sa, a = roundCorner(P1, P2, P, r)
                f_vertices, f_regions = giveCoordsCircle(O, r, sa, a, n, True)
                vertices[S1[1]] = np.array(f_vertices[0]).tolist()
                sp = (np.array(f_regions)+len(vertices)-2)[1:].tolist()
                seg += [[S1[1], sp[1][0]]]+sp[1:]
                mder = mder[1:]
                spp = copy.deepcopy(sp)
                ss1 = copy.deepcopy(S1)
                if mder:
                    mder[0][0] = spp[-1][-1]
                mascarita = mizq+[[mizq[-1][-1], ss1[1]],
                                  [ss1[1], spp[1][0]]]+spp[1:]+mder
                vertices += np.array(f_vertices)[1: -1].tolist()
                seg[fillet['end_region']][0] = len(vertices)-1
                # vertices += [O]

        original = dict(vertices=np.array(vertices), segments=np.array(seg))
        self.original = original
        if holes_dict:
            for hole in holes_dict:
                hh += [hole['center']]
                seg += (np.array(hole['regions'])+len(vertices)).tolist()
                hole['regions'] = (
                    np.array(hole['regions'])+len(vertices)).tolist()
                vertices += np.array(hole['vertices']).tolist()
            original = dict(vertices=np.array(vertices),
                            segments=np.array(seg), holes=hh)
        triangular = tr.triangulate(original, params)
        self.triangulation = triangular
        dictionary = triangular['triangles'].tolist()
        if 'o2' in params:
            tipos = ['T2V']*len(dictionary)
        else:
            tipos = ['T1V']*len(dictionary)
        gdls = triangular['vertices']
        if tipos[0] == 'T2V':
            for dicc in dictionary:
                a1 = dicc[5]
                a2 = dicc[3]
                a3 = dicc[4]
                dicc[3] = a1
                dicc[4] = a2
                dicc[5] = a3
        regions_f = []
        for s in seg:
            region = Region1D(gdls[np.ix_(s)])
            regions_f.append(region)
        Geometry2D.__init__(self, dictionary, gdls, tipos,
                            nvn=nvn, regions=regions_f, fast=fast)
        mask = []
        for region in mascarita:
            mask += [gdls[region[0]]]
        if fillets:
            self.mask = mask

        self.holes = holes_dict
        self.fillets = fillets

    @ staticmethod
    def _strdelaunay(constrained: bool = True, delaunay: bool = True, a: float = None, q: float = None, o: int = 1) -> str:
        """Create a string for the delaunay triangulation constructor

        Args:
                constrained (bool, optional): Makes the triangulation constrained. Defaults to True.
                delaunay (bool, optional): Makes all triangles delaunay. Defaults to True.
                a (float, optional): Maximum area of triange. Defaults to None.
                q (float, optional): Minimum triangle angle <=35. Defaults to None.
                o (int, optional): Order of element if 2, quadratic elements are generated. Defaults to 1.

        Returns:
                str: A string containing the input parameters for the Delaunay1V constructor
        """
        p = ''
        if o == 2:
            o = '-o2'
        else:
            o = ''
        if constrained:
            p = 'p'
        if a == None:
            a = ''
        else:
            a = 'a'+format(a)
        D = ''
        if delaunay:
            D = 'D'
        if q == None:
            q = ''
        else:
            if type(q) == int:
                if q > 35:
                    raise "No se puede crear una triangulacion con angulos menores a 35 grados"
            q = 'q'+format(q)
        return p+a+D+q+'i'+o
