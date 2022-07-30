class Solver():
    """Base Finite Element Solver.
    """

    def __init__(self, FEMObject: 'Core'):
        """Base Finite Element Solver

        Args:
            FEMObject (Core): Finite Element Problem
        """
        self.system = FEMObject
        self.type = 'Base'
        self.solutions = []
        self.solutions_info = []

    def run(self, **kargs):
        """Solves the equation system
        """

    def setSolution(self, step=-1, elements: bool = False) -> None:
        """Sets the solution to the FEM Object.

        Args:
            step (int, optional): Number of solution. Defaults to the last solution found.
            elements (bool, optional): To pass the general solution to the domain elements. Defaults to false.
        """
        self.system.solution_info = self.solutions_info[step]

        self.system.U = self.solutions[step]
        if elements:
            for e in self.system.elements:
                e.setUe(self.system.U)
