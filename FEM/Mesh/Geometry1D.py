"""General geometry 1D class.
"""


import numpy as np
from ..Utils import isBetween
import matplotlib.pyplot as plt
from ..Elements import *
from ..Elements.E1D import *
from .Geometry import *
import re


class Geometry1D(Geometry):
    """Define geometry structure

    Args:
        dictionary (list): Matrix with element definitions. Each row is an element. The gdl are defined in columns
        gdls (list): List of domain coordinates
        types (list): Types of each element
        nvn (int, optional): Nunmber of variables per node. Defaults to 1.
    """

    def __init__(self, dictionary: list, gdls: list, types: list, nvn: int = 1) -> None:
        """Define geometry structure

        Args:
            dictionary (list): Matrix with element definitions. Each row is an element. The gdl are defined in columns
            gdls (list): List of domain coordinates
            types (list): Types of each element
            nvn (int, optional): Nunmber of variables per node. Defaults to 1.
        """

        self.nvn = nvn
        self.dictionary = dictionary
        self.elements = []
        self.gdls = gdls
        self.types = types
        self.cbe = []
        self.cbn = []
        self.ngdl = int(len(self.gdls)*self.nvn)
        self.generateElements()

    @staticmethod
    def loadmsh(filename: str) -> Geometry:
        """Load geometry from previously generated MSH file

        Args:
            filename (str): Path to msh file

        Returns:
            Geometry: Output geometry
        """
        f = open(filename, 'r')
        dicc = []
        gdls = []
        types = []
        seg = []
        cbe = []
        cbn = []
        nvn = 1
        p = list(map(int, f.readline().split('\t')))
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
            cbe += [list(map(int, f.readline().split('\t')))]
        for _ in range(p[4]):
            cbn += [list(map(int, f.readline().split('\t')))]
        nvn = p[5]
        f.close()
        print('File ' + filename + ' loaded')
        o = Geometry(dicc, gdls, types, nvn, seg)
        o.cbe = cbe
        o.cbn = cbn
        return o

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

    def show(self, texto: int = 10, bolita: int = 0, figsize: list = [17, 10]) -> None:
        """Create a geometry graph

        Args:
            texto (int, optional): Text size. Defaults to 10.
            bolita (int, optional): Node size. Defaults to 0.
            figsize (list, optional): Size of figure. Defaults to [17, 10].
        """
        pass

    def saveMesh(self, ProjectName: str) -> None:
        """Saves the geometry to a MSH file with specified name

        Args:
            ProjectName (str): Project name without extension
        """
        filename = ProjectName + '.msh'
        f = open(filename, 'w')
        p = [len(self.gdls), len(self.dictionary), len(
            self.segments), len(self.cbe), len(self.cbn), self.nvn]
        f.write('\t'.join(list(map(str, p))) + '\n')
        for e in self.gdls:
            f.write('\t'.join(list(map(str, e))) + '\n')
        for e in self.types:
            f.write(e + '\n')
        for e in self.dictionary:
            f.write('\t'.join(list(map(str, e))) + '\n')
        for e in self.segments:
            f.write('\t'.join(list(map(str, e))) + '\n')
        for e in self.cbe:
            f.write('\t'.join(list(map(str, e))) + '\n')
        for e in self.cbn:
            f.write('\t'.join(list(map(str, e))) + '\n')
        f.close()
        print('File ' + filename + ' saved')
