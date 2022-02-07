"""Defines a Lagrangian 2 order triangular element 
"""


from .Element2D import Element2D, np
from .TriangularScheme import TriangularScheme
from ..E1D.QuadraticElement import QuadraticElement


class QTriangular(Element2D, TriangularScheme):
    """Creates a lagrangian element of order 2

    Args:
        coords (np.ndarray): Element coordinates matrix
        gdl (np.ndarray): Element gdl matrix
        n (int, optional): Number of Gauss Points. Defaults to 2.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 3, **kargs) -> None:
        """Creates a lagrangian element of order 2

        Args:
            coords (np.ndarray): Element coordinates matrix
            gdl (np.ndarray): Element gdl matrix
            n (int, optional): Number of Gauss Points. Defaults to 2.
        """

        coords = np.array(coords)
        delta = coords[0]-coords[1]
        he1 = np.linalg.norm(delta)
        e1 = QuadraticElement(np.array([[0], [he1*0.5], [he1]]),
                              np.array([[-1, -1, -1]]), border=True)

        delta = coords[2]-coords[1]
        he2 = np.linalg.norm(delta)
        e2 = QuadraticElement(np.array([[0], [he2*0.5], [he2]]),
                              np.array([[-1, -1, -1]]), border=True)

        delta = coords[0]-coords[2]
        he3 = np.linalg.norm(coords[0]-coords[2])
        e3 = QuadraticElement(np.array([[0], [he3*0.5], [he3]]),
                              np.array([[-1, -1, -1]]), border=True)
        self.borders = [e1, e2, e3]

        _coords = np.array([coords[i] for i in range(3)])
        Element2D.__init__(self, coords, _coords, gdl, **kargs)
        TriangularScheme.__init__(self, n)

    @classmethod
    def description(self):
        self = TriangularScheme(3)

        def bmatrix(a, header=[], caption='', table=False):
            rv = []
            if len(a.shape) > 2:
                raise ValueError('bmatrix can at most display two dimensions')
            if table:
                rv += [r"""\begin{table}
        \centering
        \caption{"""+caption+"""}"""]
                rv += [r'\begin{tabular}{' + '|c'*len(a[0]) + '|}\hline']
                rv += ['  ' +
                       ' & '.join([r'\textbf{'+i+'}' for i in header])+r'\\\hline']
                lines = str(a).replace('[', '').replace(
                    ']', '').splitlines()
                rv += ['  ' + ' & '.join(l.split()) +
                       r'\\ \hline' for l in lines]
                rv += [r'\end{tabular}']
                rv += [r"""\end{table}"""]
            else:
                lines = str(a).replace('[', '').replace(
                    ']', '').splitlines()
                rv += [r'\begin{bmatrix}']
                rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
                rv += [r'\end{bmatrix}']
            return '\n'.join(rv)
        return r"""Elemento de 6 nodos. El primer nodo se encuentra centrado en el origen y los otros dos nodos principales se encuentran a una unidad en ambas direcciones. Se crean nodos secundarios en los segmentos que forma cada pareja de nodos primaros. Los nodos se enumeran en el orden contrario a las manecillas del reloj empezando por los nodos primarios y posteriormente los nodos secundarios.
        Para calcular las integrales de este elemento se usan """+format(len(self.Z))+r""" puntos de Gauss:\\$$\zeta="""+bmatrix(self.Z)+"""$$$$W="""+bmatrix(self.W)+r"""$$\\\\ Para este elemento se usaron las siguientes funciones de forma:
        $$\Psi_0=2(\zeta+\eta-1)(\zeta+\eta-0.5)$$
        $$\Psi_1=2\zeta(\zeta-0.5)$$
        $$\Psi_2=2\eta(\eta-0.5)$$
        $$\Psi_3=-4(\zeta+\eta-1)(\zeta)$$
        $$\Psi_4=4\zeta\eta$$
        $$\Psi_5=-4\eta(\zeta+\eta-1)$$
        """

    def psis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """
        return np.array([
            2.0*(z[0]+z[1]-1.0)*(z[0]+z[1]-0.5),
            2.0*z[0]*(z[0]-0.5),
            2.0*z[1]*(z[1]-0.5),
            -4.0*(z[0]+z[1]-1.0)*(z[0]),
            4.0*z[0]*z[1],
            -4.0*z[1]*(z[0]+z[1]-1.0)]).T

    def dpsis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """

        return np.array([
            [4.0*z[0]+4.0*z[1]-3.0, 4.0*z[1]+4.0*z[0]-3.0],
            [4.0*z[0]-1.0, 0*z[0]],
            [0*z[0], 4.0*z[1]-1.0],
            [-8.0*z[0]-4.0*(z[1]-1.0), -4.0*z[0]],
            [4.0*z[1], 4.0*z[0]],
            [-4.0*z[1], -8.0*z[1]-4.0*z[0]+4.0]
        ])
