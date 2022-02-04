"""Defines a Lagrangian 1 order triangular element 
"""


from .Element2D import Element2D, np
from ..E1D.LinealElement import LinealElement
from .TriangularScheme import TriangularScheme


class LTriangular(Element2D, TriangularScheme):
    """Creates a lagrangian triangular element of order 1

    Args:
        coords (np.ndarray): Element coordinates matrix
        gdl (np.ndarray): Element gdl matrix
        n (int, optional): Number of Gauss Points. Defaults to 2.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 2, **kargs) -> None:
        """Creates a lagrangian triangular element of order 1

        Args:
            coords (np.ndarray): Element coordinates matrix
            gdl (np.ndarray): Element gdl matrix
            n (int, optional): Number of Gauss Points. Defaults to 2.
        """

        coords = np.array(coords)
        he1 = np.linalg.norm(coords[1]-coords[0])
        e1 = LinealElement(np.array([[0], [he1]]),
                           np.array([[-1, -1]]), border=True)

        he2 = np.linalg.norm(coords[2]-coords[1])
        e2 = LinealElement(np.array([[0], [he2]]),
                           np.array([[-1, -1]]), border=True)

        he3 = np.linalg.norm(coords[0]-coords[2])
        e3 = LinealElement(np.array([[0], [he3]]),
                           np.array([[-1, -1]]), border=True)

        self.borders = [e1, e2, e3]
        Element2D.__init__(self, coords, coords, gdl, **kargs)
        TriangularScheme.__init__(self, n)

    @classmethod
    def description(self):
        self = TriangularScheme(2)

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
        return r"""Elemento de 3 nodos. El primer nodo se encuentra centrado en el origen y los otros dos nodos se encuentran a una unidad en ambas direcciones. Los nodos se enumeran en el orden contrario a las manecillas del reloj.
        Para calcular las integrales de este elemento se usan """+format(len(self.Z))+r""" puntos de Gauss:\\$$\zeta="""+bmatrix(self.Z)+"""$$$$W="""+bmatrix(self.W)+r"""$$\\\\ Para este elemento se usaron las siguientes funciones de forma:
        $$\Psi_0=1-\zeta-\eta$$
        $$\Psi_1=\zeta$$
        $$\Psi_2=\eta$$
        """

    def psis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """

        return np.array([
            1.0-z[0]-z[1],
            z[0],
            z[1]]).T

    def dpsis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """

        kernell = (z[0]-z[0])
        return np.array([
            [-1.0*(1+kernell), -1.0*(1+kernell)],
            [1.0*(1+kernell), 0.0*(1+kernell)],
            [0.0*(1+kernell), 1.0*(1+kernell)]
        ])
