"""Define the structure of a finite element problem. This is the parent class of individual equation classes.

The individual children classes must implement the method for calculating the element matrices and post processing.
"""


from typing import Union, Callable
from tqdm import tqdm
import numpy as np
from .Geometry import Geometry
from .Solvers import Lineal, NonLinealSolver, LinealSparse, Solver, Parabolic
import logging
from .FEMLogger import FEMLogger
from functools import partialmethod
from .Elements import Element
import json


class Core():
    """Create the Finite Element problem.

        Args:
            geometry (Geometry): Input geometry. The geometry must contain the elements, and the boundary conditions.
            You can create the geometry of the problem using the Geometry class.
            solver (Union[Lineal, NonLinealSolver], optional): Finite Element solver. If not provided, Lineal solver is used.
            sparse (bool, optional): To use sparse matrix formulation. Defaults to False
            verbose (bool, optional): To print console messages and progress bars. Defaults to False.

    """

    def __init__(self, geometry: Geometry, solver: Union[Lineal, NonLinealSolver] = None, sparse: bool = False, verbose: bool = False, name='') -> None:
        """Create the Finite Element problem.

            Args:
                geometry (Geometry): Input geometry. The geometry must contain the elements, and the boundary conditions.
                You can create the geometry of the problem using the Geometry class.
                solver (Union[Lineal, NonLinealSolver], optional): Finite Element solver. If not provided, Lineal solver is used.
                sparse (bool, optional): To use sparse matrix formulation. Defaults to False
                verbose (bool, optional): To print console messages and progress bars. Defaults to False.
                name (str, optional): To print custom name on logging file. Defaults to ''.

        """
        self.logger: FEMLogger = FEMLogger(name)
        if verbose:
            self.logger.setup_logging(console_log_level="info")
        else:
            self.logger.setup_logging()
        self.geometry: Geometry = geometry
        self.ngdl: int = self.geometry.ngdl
        if not sparse:
            self.K: np.ndarray = np.zeros([self.ngdl, self.ngdl])
        self.F: np.ndarray = np.zeros([self.ngdl, 1])
        self.Q: np.ndarray = np.zeros([self.ngdl, 1])
        self.U: np.ndarray = np.zeros([self.ngdl, 1])
        self.S: np.ndarray = np.zeros([self.ngdl, 1])
        self.cbe: list = self.geometry.cbe
        self.cbn: list = self.geometry.cbn
        self.elements: list[Element] = self.geometry.elements
        self.verbose: bool = verbose
        tqdm.__init__ = partialmethod(tqdm.__init__, disable=not verbose)
        self.name: str = 'Generic FEM '

        if not solver:
            self.solver: Solver = Lineal(self)
            if sparse:
                self.solver: Solver = LinealSparse(self)
        else:
            self.solver: Solver = solver(self)
        if self.solver.type == 'non-lineal-newton':
            self.T: np.ndarray = np.zeros([self.ngdl, self.ngdl])
        elif self.solver.type == "Base":
            logging.error("Base solver should not be used.")
            raise Exception("Base solver should not be used.")
        self.properties: dict = {'verbose': self.verbose,
                                 'name': name, 'problem': self.__class__.__name__}

    def description(self):
        """Generates the problem description for loggin porpuses
        """
        return f'FEM problem using the {self.name} formulation.,DOF: {self.ngdl},Elements: {len(self.elements)},Solver: {self.solver.type},EBC: {len(self.cbe)},NBC: {len(self.cbn)}'

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

    def restartMatrix(self) -> None:
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

    def boundaryConditions(self) -> None:
        """Assign boundary conditions to the system. 
        The boundary conditios are assigned in this order:

        1. Natural boundary conditions
        2. Essential boundary conditions

        This ensures that in a node with 2 boundary conditions
        the essential boundary conditions will be applied.
        """
        logging.info('Boundary conditions...')
        for i in tqdm(self.cbn, unit=' Natural'):
            self.Q[int(i[0])] = i[1]

        if self.cbe:
            boundary_conditions = np.zeros([self.ngdl, 1])
            cb = np.array(self.cbe)
            ncb = len(cb)
            boundary_conditions[np.ix_(cb[:, 0].astype(int))
                                ] = cb[:, 1].reshape([ncb, 1])
            # FIXME esto puede resultar en pasar la matriz dwe lil a densa!!
            self.S = self.S - (boundary_conditions.T@self.K).T
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

    def condensedSystem(self) -> None:
        """Assign boundary conditions to the system and modifies the matrices to condensed mode 
        The boundary conditios are assigned in this order:

        1. Natural boundary conditions
        2. Essential boundary conditions

        This ensures that in a node with 2 boundary conditions
        the essential boundary conditions will be applied.
        """
        # FIXME tiene que funcionar para todos los casos
        # self.boundaryConditions()
        logging.info('Boundary conditions...')
        for i in tqdm(self.cbn, unit=' Natural'):
            self.Q[int(i[0])] = i[1]

        if self.cbe:
            boundary_conditions = np.zeros([self.ngdl, 1])
            cb = np.array(self.cbe)
            ncb = len(cb)
            boundary_conditions[np.ix_(cb[:, 0].astype(int))
                                ] = cb[:, 1].reshape([ncb, 1])
            # FIXME esto puede resultar en pasar la matriz dwe lil a densa!!
            self.S = self.S - (boundary_conditions.T@self.K).T
            for i in tqdm(self.cbe, unit=' Essential'):
                self.K[int(i[0]), :] = 0
                self.K[:, int(i[0])] = 0
                self.K[int(i[0]), int(i[0])] = 1
                if self.calculateMass:
                    self.M[int(i[0]), :] = 0
                    self.M[:, int(i[0])] = 0
                    self.M[int(i[0]), int(i[0])] = 1
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
        desc = self.description().split(',')
        for d in desc:
            logging.debug(d)
        self.solveES(**kargs)
        if plot:
            logging.info('Post processing solution...')
            self.postProcess(**kargs)
            logging.info('Done!')
        duration = self.logger.end_timer().total_seconds()
        self.properties['duration'] = duration
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

    def exportJSON(self, filename: str = None):
        self.properties['description'] = self.description()
        if not filename:
            filename = __name__
        y = self.geometry.exportJSON()
        pjson = json.loads(y)
        pjson["solutions"] = []
        # Avoid Numpy array conversion
        for info, sol in zip(self.solver.solutions_info, np.array(self.solver.solutions).tolist()):
            pjson["solutions"].append({'info': info, "U": sol})
        pjson['properties'] = self.properties
        y = json.dumps(pjson)
        with open(filename, "w") as f:
            f.write(y)


class CoreTransient(Core):
    """docstring for CoreTransient
    """

    def __init__(self, geometry: Geometry, solver: Solver = None, verbose: bool = False, name='') -> None:
        Core.__init__(self, geometry=geometry, solver=solver,
                      verbose=verbose, name=name, sparse=False)
        self.dt: float = 0.1
        self.t: float = 0.0

        # Initial condition for all gdls in domain
        self.u0: list[float] = [0.0]*self.ngdl

    def set_initial_condition(self, ic: Union[list[float], Callable[[], float], float]) -> None:
        if isinstance(ic, float):
            self.u0 = [ic]*self.ngdl
        elif isinstance(ic, Callable):
            self.u0 = ic(self.geometry.gdls)
        elif isinstance(ic, list):
            self.u0 = ic

    def apply_initial_condition(self, t0, dt):
        self.U[:, 0] = self.u0
        for i in self.cbe:
            self.U[int(i[0]), 0] = i[1]
        self.solver.solutions.append(self.U.copy())
        self.solver.solutions_info.append(
            {'solver-type': self.solver.type, "time": t0, "dt": dt})
        self.solver.setSolution(-1, True)

    def ensembling(self) -> None:
        logging.info('Ensembling equation system...')
        for e in self.elements:
            self.K[np.ix_(e.gdlm, e.gdlm)] += e.Ke
            self.F[np.ix_(e.gdlm)] += e.Fe
        logging.info('Done!')


class CoreParabolic(CoreTransient):
    """docstring for CoreParabolic
    """

    def __init__(self, geometry: Geometry, solver: Solver = None, verbose: bool = False, name=''):
        if not solver:
            solver = Parabolic
        CoreTransient.__init__(self, geometry=geometry, solver=solver,
                               verbose=verbose, name=name)

        self.alpha: float = 0.5  # Crack nocholson. Subclases maybe???

    def set_alpha(self, alpha):
        self.alpha = alpha
        # TODO select solver acording to alpha value

    def postProcess(self, t0, tf, steps) -> None:
        """Post process the solution
        """
        pass


class CoreHiperbolic(CoreTransient):
    """docstring for CoreParabolic
    """

    def __init__(self, geometry: Geometry, solver: Solver = None, verbose: bool = False, name=''):
        CoreTransient.__init__(self, geometry=geometry, solver=solver,
                               verbose=verbose, name=name)
        self.U_dot: np.ndarray = np.zeros([self.ngdl, 1])
        self.U_dot_dot: np.ndarray = np.zeros([self.ngdl, 1])
        self.du0: list[float] = [0.0]*self.ngdl
        self.ddu0: list[float] = [0.0]*self.ngdl
        self.dU: np.ndarray = np.zeros([self.ngdl, 1])
        self.ddU: np.ndarray = np.zeros([self.ngdl, 1])
        self.alpha: float = 0.5
        self.gamma: float = 0.5
        self.apply_initial_condition()

    def set_initial_condition(self, u0: Union[list[float], Callable[[], float], float], du0: Union[list[float], Callable[[], float], float], ddu0: Union[list[float], Callable[[], float], float]) -> None:
        if isinstance(u0, float):
            self.u0 = [u0]*self.ngdl
        elif isinstance(u0, Callable):
            self.u0 = u0(self.geometry.gdls)
        elif isinstance(u0, list):
            self.u0 = u0

        if isinstance(du0, float):
            self.du0 = [du0]*self.ngdl
        elif isinstance(du0, Callable):
            self.du0 = du0(self.geometry.gdls)
        elif isinstance(du0, list):
            self.du0 = du0

        if isinstance(ddu0, float):
            self.ddu0 = [ddu0]*self.ngdl
        elif isinstance(ddu0, Callable):
            self.ddu0 = ddu0(self.geometry.gdls)
        elif isinstance(ddu0, list):
            self.ddu0 = ddu0

    def apply_initial_condition(self):
        self.U[:, 0] = self.u0
        self.dU[:, 0] = self.du0
        self.ddU[:, 0] = self.ddu0
        self.solver.setSolution(-1, True)
