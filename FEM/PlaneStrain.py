"""2D Elasticity: Plane Strain
"""


from .PlaneStress import PlaneStress, Geometry
from typing import Tuple, Callable


class PlaneStrain(PlaneStress):
    """Create a Plain Strain problem

    Args:
            geometry (Geometry): 2D 2 variables per node geometry
            E (int or float or list): Young Moduli. If number, all element will have the same young moduli. If list, each position will be the element young moduli, so len(E) == len(self.elements)
            v (int or float or list): Poisson ratio. If number, all element will have the same Poisson ratio. If list, each position will be the element Poisson ratio, so len(v) == len(self.elements)
            fx (function, optional): Function fx, if fx is constant you can use fx = lambda x: [value]. Defaults to lambda x:0.
            fy (function, optional): Function fy, if fy is constant you can use fy = lambda x: [value]. Defaults to lambda x:0.
    """

    def __init__(self, geometry: Geometry, E: Tuple[float, list], v: Tuple[float, list], fx: Callable = lambda x: 0, fy: Callable = lambda x: 0) -> None:
        """Create a Plain Strain problem

        Args:
                geometry (Geometry): 2D 2 variables per node geometry
                E (int or float or list): Young Moduli. If number, all element will have the same young moduli. If list, each position will be the element young moduli, so len(E) == len(self.elements)
                v (int or float or list): Poisson ratio. If number, all element will have the same Poisson ratio. If list, each position will be the element Poisson ratio, so len(v) == len(self.elements)
                fx (function, optional): Function fx, if fx is constant you can use fx = lambda x: [value]. Defaults to lambda x:0.
                fy (function, optional): Function fy, if fy is constant you can use fy = lambda x: [value]. Defaults to lambda x:0.
        """

        PlaneStress.__init__(self, geometry, E, v, 1, fx, fy)
        self.C11 = []
        self.C12 = []
        self.C66 = []
        for i in range(len(self.E)):
            C11 = self.E[i]*(1-self.v[i])/(1+self.v[i])/(1-2*self.v[i])
            C12 = self.E[i]*(self.v[i])/(1+self.v[i])/(1-2*self.v[i])
            C66 = self.E[i] / 2 / (1 + self.v[i])
            self.C11.append(C11)
            self.C12.append(C12)
            self.C66.append(C66)
