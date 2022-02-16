"""General geometry class.
"""

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
from typing import Callable
from ast import literal_eval, parse
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

    def __init__(self, dictionary: list, gdls: list, types: list, nvn: int = 1, segments: list = [], fast=False) -> None:
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
        self.nvn = nvn
        self.dictionary = dictionary
        self.elements = []
        self.gdls = gdls
        self.types = types
        self.segments = segments
        self.cbe = []
        self.cbn = []
        self.centroids = []
        self.fast = fast
        self.initialize()
        self.calculateCentroids()

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
        print('Detecting non local elements')
        #TODO hacerlo mas eficiente.
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

    @staticmethod
    def loadmsh(filename: str, **kargs):
        """Load geometry from previously generated MSH file

        Args:
            filename (str): Path to msh file

        Returns:
            Geometry: Output geometry
        """
        #TODO poner esto como un constructor?
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
        if len_mask == 0:
            mask = None
        for _ in range(len_holes):
            holes.append(literal_eval(f.readline()))
        for _ in range(len_fillets):
            fillets.append(literal_eval(f.readline()))
        for _ in range(len_mask):
            mask += [list(map(float, f.readline().split('\t')))]
        f.close()
        print('File ' + filename + ' loaded')
        o = Geometry(dicc, gdls, types, nvn, seg, **kargs)
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
        self.elements = [0.0]*len(self.dictionary)
        for i, d in enumerate(tqdm(self.dictionary, unit='Element')):
            coords = np.array(self.gdls)[np.ix_(d)]
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

            # if i > 1 and isinstance(element, self.elements[-1].__class__):
            #     _p = self.elements[-1]._p
            #     dpz = self.elements[-1].dpz
            # else:
            #     _p = element.psis(element.Z.T)
            #     dpz = element.dpsis(element.Z.T).T
            # element.fastInit(_p, dpz)

            self.elements[i] = element
        print('Done!')

    def show(self):
        """Creates a geometry graph"""
        pass

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

    def calculateCentroids(self) -> None:
        """Calculate elements centroids
        """
        for e in self.elements:
            x,p=e.T(e.center.T)
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

    def giveNodesOfSegment(self, segment: int, tol: float) -> np.ndarray:
        """Give nodes over a segment

        Args:
            segment (int): Segment number. Start with 0
            tol (float): Tolerance for finding nearest nodes

        Returns:
            np.ndarray: List of nodes in the specified segment
        """
        #TODO Hacer que esto sea n dimensional
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
            list: List of elements in the specified segment
        """
        a = []
        nodes = self.giveNodesOfSegment(segment, tol)
        for e in self.elements:
            if np.sum(np.isin(e.gdl[0], nodes*self.nvn)) > 0:
                a.append(e)
        return a

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

    def cbeAllBorders(self, value: float, tol: float = 1*10**(-5)) -> None:
        """Set all segments border conditions to the specified value

        Args:
            value (float): Value of the border condition
            tol (float, optional): Tolerance for finding near nodes in segments. Defaults to 1*10**(-5).
        """
        for s in range(len(self.segments)):
            for i in range(self.nvn):
                self.cbe += self.cbFromSegment(s, value, (i+1), tol)

    def exportJSON(self, filename: str = None) -> str:
        """Export geometry definition as JSON file or JSON string

        Args:
            filename (str, optional): If given, a JSON file is created. Defaults to None.

        Returns:
            str: JSON string
        """
        #TODO Faltan los holes y los fillets
        x = {
            "nodes": self.gdls,
            "dictionary": self.dictionary,
            "types": self.types,
            "regions": self.segments,
            "ebc": self.cbe,
            "nbc": self.cbn,
            "nvn": self.nvn,
            "ngdl": self.ngdl,
        }
        y = json.dumps(x)
        if filename:
            with open(filename, "w") as f:
                f.write(y)
        return y

    @staticmethod
    def importJSON(filename: str, **kargs) -> 'Geometry':
        """Import geometry definition from JSON file

        Args:
            filename (str): Path to the JSON file

        Returns:
            Geometry: Generated JSON file
        """
        with open(filename) as f:
            parsed = json.loads(f.read())
            dcc = parsed['dictionary']
            nodes = parsed['nodes']
            types = parsed['types']
            nvn = parsed['nvn']
            regions = parsed['regions']
            o = Geometry(dcc, nodes, types, nvn, regions, **kargs)
            o.cbe = parsed['ebc']
            o.cbn = parsed['nbc']
            return o
