"""Define the structure of a finite element problem. This is the parent class of individual equation classes.

The individual children classes must implement the method for calculating the element matrices and post processing.
"""


from asyncio.log import logger
from typing import Union
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from .Geometry import Geometry
import os
from .Elements import *
from .Solvers import *
import logging
from .FEMLogger import FEMLogger
from functools import partialmethod


class Core():
    def __init__(self, geometry: Geometry, solver: Union[Lineal, NonLinealSolver] = None, sparse: bool = False, verbose: bool = False) -> None:
        """Create the Finite Element problem.

            Args:
                geometry (Geometry): Input geometry. The geometry must contain the elements, and the border conditions.
                You can create the geometry of the problem using the Geometry class.
                solver (Union[Lineal, NonLinealSolver], optional): Finite Element solver. If not provided, Lineal solver is used.
                sparse (bool, optional): To use sparse matrix formulation. Defaults to False
                verbose (bool, optional): To print console messages and progress bars. Defaults to False.

        """
        self.logger = FEMLogger(__name__)
        if verbose:
            self.logger.setup_logging(console_log_level="info")
        else:
            self.logger.setup_logging()
        self.geometry = geometry
        self.ngdl = self.geometry.ngdl
        if not sparse:
            self.K = np.zeros([self.ngdl, self.ngdl])
        self.F = np.zeros([self.ngdl, 1])
        self.Q = np.zeros([self.ngdl, 1])
        self.U = np.zeros([self.ngdl, 1])
        self.S = np.zeros([self.ngdl, 1])
        self.cbe = self.geometry.cbe
        self.cbn = self.geometry.cbn
        self.elements = self.geometry.elements
        self.verbose = verbose
        tqdm.__init__ = partialmethod(tqdm.__init__, disable=not verbose)

        if not solver:
            self.solver = Lineal(self)
        else:
            self.solver = solver(self)
        if self.solver.type == 'non-lineal-newton':
            self.T = np.zeros([self.ngdl, self.ngdl])

    def ensembling(self) -> None:
        """Ensembling of equation system. This method use the element gdl
        and the element matrices. The element matrices degrees of fredom must
        match the dimension of the element gdl. For m>1 variables per node,
        the gdl will be flattened. This ensure that the element matrices will always
        be a 2-D Numpy Array.
        """
        logging.info('Ensembling equation system...')

        for e in self.elements:
            self.K[np.ix_(e.gdlm, e.gdlm)] += e.Ke
            self.F[np.ix_(e.gdlm)] += e.Fe
            self.Q[np.ix_(e.gdlm)] += e.Qe
            if 'newton' in self.solver.type:
                try:
                    self.T[np.ix_(e.gdlm, e.gdlm)] += e.Te
                except Exception as e:
                    logging.error(
                        "Impossible to access tangent matrix. Check tangent matrix creation in integration class.")
                    raise e
        logging.info('Done!')

    def restartMatrix(self):
        """Sets all model matrices and vectors to 0 state
        """
        self.K[:, :] = 0.0
        self.F[:, :] = 0.0
        self.Q[:, :] = 0.0
        self.S[:, :] = 0.0
        if 'newton' in self.solver.type:
            try:
                self.T[:, :] = 0.0
            except Exception as e:
                logging.error("Impossible to clear tangent matrix.")
                raise e

    def borderConditions(self) -> None:
        """Assign border conditions to the system. 
        The border conditios are assigned in this order:

        1. Natural border conditions
        2. Essential border conditions

        This ensures that in a node with 2 border conditions
        the essential border conditions will be applied.
        """
        logging.info('Border conditions...')
        for i in tqdm(self.cbn, unit=' Natural'):
            self.Q[int(i[0])] = i[1]

        if self.cbe:
            border_conditions = np.zeros([self.ngdl, 1])
            cb = np.array(self.cbe)
            ncb = len(cb)
            border_conditions[np.ix_(cb[:, 0].astype(int))
                              ] = cb[:, 1].reshape([ncb, 1])
            self.S = self.S - (border_conditions.T@self.K).T
            for i in tqdm(self.cbe, unit=' Essential'):
                self.K[int(i[0]), :] = 0
                self.K[:, int(i[0])] = 0
                self.K[int(i[0]), int(i[0])] = 1
                if 'newton' in self.solver.type:
                    try:
                        self.T[int(i[0]), :] = 0
                        self.T[:, int(i[0])] = 0
                        self.T[int(i[0]), int(i[0])] = 1
                    except Exception as e:
                        logging.error("Impossible to access tangent matrix.")
                        raise e

        self.S = self.S + self.F + self.Q
        for i in self.cbe:
            self.S[int(i[0])] = i[1]
        logging.info('Done!')

    def solveES(self, **kargs) -> None:
        """Solve the finite element problem
        """
        self.solver.run(**kargs)

    def solve(self, plot: bool = True, **kargs) -> None:
        """A series of Finite Element steps

        Args:
            plot (bool, optional): To post process the solution. Defaults to True.
            **kargs: Solver specific parameters.
        """
        if self.solver.type == 'lineal':
            logging.info('Creating element matrices...')
            self.elementMatrices()
            logging.info('Done!')
            self.ensembling()
            self.borderConditions()
        self.solveES(**kargs)
        if plot:
            logging.info('Post processing solution...')
            self.postProcess(**kargs)
            logging.info('Done!')
        self.logger.end_timer()
        logging.info("End!")

    def solveFromFile(self, file: str, plot: bool = True, **kargs) -> None:
        """Load a solution file and show the post process for a given geometry

        Args:
                file (str): Path to the previously generated solution file.
                plot (bool, optional): To post process the solution. Defaults to True.
        """
        logging.info('Loading File...')
        self.U = np.loadtxt(file)
        for e in self.elements:
            e.setUe(self.U)
        logging.info('Done!')
        if plot:
            logging.info('Post processing solution...')
            self.postProcess(**kargs)
            logging.info('Done!')

    def solveFromArray(self, solution: np.ndarray, plot: bool = True, **kargs) -> None:
        """Load a solution array to the problem.

        Args:
                solution (np.ndarray): Solution vertical array with shape (self.ngdl,1)
                plot (bool, optional): To post process the solution. Defaults to True.
        """
        logging.info('Casting solution')
        self.U = solution
        for e in self.elements:
            e.setUe(self.U)
        logging.info('Done!')
        if plot:
            logging.info('Post processing solution...')
            self.postProcess(**kargs)
            logging.info('Done!')

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
