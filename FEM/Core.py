"""Define the structure of a finite element problem. This is the parent class of individual equation classes.

The individual children classes must implement the method for calculating the element matrices and post processing.
"""


from tqdm import tqdm
import numpy as np
from .Mesh import Geometry


class Core():
    def __init__(self, geometry: Geometry) -> None:
        """Create the Finite Element problem.

            Args:
                geometry (Geometry): Input geometry. The geometry must contain the elements, and the border conditions.
                You can create the geometry of the problem using the Geometry class.
        """
        self.geometry = geometry
        self.ngdl = self.geometry.ngdl
        self.K = np.zeros([self.ngdl, self.ngdl])
        self.F = np.zeros([self.ngdl, 1])
        self.Q = np.zeros([self.ngdl, 1])
        self.U = np.zeros([self.ngdl, 1])
        self.S = np.zeros([self.ngdl, 1])
        self.cbe = self.geometry.cbe
        self.cbn = self.geometry.cbn
        self.elements = self.geometry.elements

    def ensembling(self) -> None:
        """Ensembling of equation system. This method use the element gdl
        and the element matrices. The element matrices degrees of fredom must
        match the dimension of the element gdl. For m>1 variables per node,
        the gdl will be flattened. This ensure that the element matrices will always 
        be a 2-D Numpy Array.
        """
        print('Ensembling equation system...')
        for e in tqdm(self.elements, unit='Element'):
            self.K[np.ix_(e.gdlm, e.gdlm)] += e.Ke
            self.F[np.ix_(e.gdlm)] += e.Fe
            self.Q[np.ix_(e.gdlm)] += e.Qe
        self._K = np.copy(self.K)
        print('Done!')

    def borderConditions(self) -> None:
        """Assign border conditions to the system. 
        The border conditios are assigned in this order:

        1. Natural border conditions
        2. Essential border conditions

        This ensures that in a node with 2 border conditions
        the essential border conditions will be applied.
        """
        print('Border conditions...')
        for i in tqdm(self.cbn, unit=' Natural'):
            self.Q[int(i[0])] = i[1]
        for i in tqdm(self.cbe, unit=' Essential'):
            ui = np.zeros([self.ngdl, 1])
            ui[int(i[0])] = i[1]
            vv = np.dot(self.K, ui)
            self.S = self.S - vv
            self.K[int(i[0]), :] = 0
            self.K[:, int(i[0])] = 0
            self.K[int(i[0]), int(i[0])] = 1
        self.S = self.S + self.F + self.Q
        for i in self.cbe:
            self.S[int(i[0])] = i[1]
        print('Done!')

    def solveES(self, path: str = '') -> None:
        """Solve the equation system using numpy.solve algorithm

        Args:
                path (str, optional): Path to save a text file with the solution of the problem
                This file can be loaded witouth spendign time in other finite element steps. Defaults to ''.
        """
        print('Solving equation system...')
        self.U = np.linalg.solve(self.K, self.S)
        if not path == '':
            np.savetxt(path, self.U, delimiter=',')
        for e in self.elements:
            e.setUe(self.U)
        print('Done!')

    def solve(self, path: str = '', plot: bool = True) -> None:
        """A series of Finite Element steps

        Args:
                path (str, optional): Path to save a text file with the solution of the problem
                This file can be loaded witouth spendign time in other finite element steps. Defaults to ''.
        """
        print('Creating element matrices...')
        self.elementMatrices()
        print('Done!')
        self.ensembling()
        self.borderConditions()
        self.solveES(path)
        if plot:
            print('Post processing solution...')
            self.postProcess()
            print('Done!')

    def solveFromFile(self, file: str, plot: bool = True) -> None:
        """Load a solution file and show the post process for a given geometry

        Args:
                file (str): Path to the previously generated solution file.
        """
        print('Loading File...')
        self.U = np.loadtxt(file)
        for e in self.elements:
            e.setUe(self.U)
        print('Done!')
        if plot:
            print('Post processing solution...')
            self.postProcess()
            print('Done!')

    def profile(self) -> None:
        """Create a profile for a 3D or 2D problem.
        """
        pass

    def elementMatrices(self) -> None:
        """Calculate the element matrices
        """
        pass

    def postProcess(self) -> None:
        """Post process the solution
        """
        pass
