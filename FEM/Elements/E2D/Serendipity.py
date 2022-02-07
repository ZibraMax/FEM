"""Define a o2 serendipity element"""


from .Element2D import Element2D, np
from .RectangularScheme import RectangularScheme
from ..E1D.QuadraticElement import QuadraticElement


class Serendipity(Element2D, RectangularScheme):
    """Creates a Serendipity element

    Args:
        coords (np.ndarray): Coordinate matrix of element
        gdl (np.ndarray): Coordinate matrix of element GDL's
        n (int, optional): Number of gauss points. Defaults to 3.
    """

    def __init__(self, coords: np.ndarray, gdl: np.ndarray, n: int = 3, **kargs) -> None:
        """Creates a Serendipity element

        Args:
            coords (np.ndarray): Coordinate matrix of element
            gdl (np.ndarray): Coordinate matrix of element GDL's
            n (int, optional): Number of gauss points. Defaults to 3.
        """

        _coords = np.array([coords[i] for i in range(4)])
        coords = np.array(coords)

        he1 = np.linalg.norm(coords[1]-coords[0])
        e1 = QuadraticElement(np.array([[0], [he1*0.5], [he1]]),
                              np.array([[-1, -1, -1]]), border=True)

        he2 = np.linalg.norm(coords[2]-coords[1])
        e2 = QuadraticElement(np.array([[0], [he2*0.5], [he2]]),
                              np.array([[-1, -1, -1]]), border=True)

        he3 = np.linalg.norm(coords[3]-coords[2])
        e3 = QuadraticElement(np.array([[0], [he3*0.5], [he3]]),
                              np.array([[-1, -1, -1]]), border=True)

        he4 = np.linalg.norm(coords[0]-coords[3])
        e4 = QuadraticElement(np.array([[0], [he4*0.5], [he4]]),
                              np.array([[-1, -1, -1]]), border=True)

        self.borders = [e1, e2, e3, e4]
        Element2D.__init__(self, coords, _coords, gdl, **kargs)
        RectangularScheme.__init__(self, n)

    @classmethod
    def description(self):
        self = RectangularScheme(3)

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
        return r"""Elemento de 8 nodos. El primer nodo se encuentra en las coordenadas $\zeta,\eta=[-1,-1]$ y los otros nodos principales se encuentran a dos unidadades en ambas direcciones. Se crean nodos secundarios en los segmentos que forma cada pareja de nodos primaros. Los nodos se enumeran en el orden contrario a las manecillas del reloj empezando por los nodos primarios y posteriormente los nodos secundarios.
        Para calcular las integrales de este elemento se usan """+format(len(self.Z))+r""" puntos de Gauss:\\$$\zeta="""+bmatrix(self.Z)+"""$$$$W="""+bmatrix(self.W)+r"""$$\\\\ Para este elemento se usaron las siguientes funciones de forma:
        $$\Psi_0=0.25(1-\zeta)(1-\eta)(-1-\zeta-\eta)$$
        $$\Psi_1=0.25(1+\zeta)(1-\eta)(-1+\zeta-\eta)$$
        $$\Psi_2=0.25(1+\zeta)(1+\eta)(-1+\zeta+\eta)$$
        $$\Psi_3=0.25(1-\zeta)*(1+\eta)*(-1-\zeta+\eta)$$
        $$\Psi_4=0.5(1-\zeta^2)(1-\eta)$$
        $$\Psi_5=0.5(1+\zeta)(1-\eta^2)$$
        $$\Psi_6=0.5(1-\zeta^2)(1+\eta)$$
        $$\Psi_7=0.5(1-\zeta)(1-\eta^2)$$
        """

    def psis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function evaluated in Z points
        """

        return np.array([
            0.25*(1.0-z[0])*(1.0-z[1])*(-1.0-z[0]-z[1]),
            0.25*(1.0+z[0])*(1.0-z[1])*(-1.0+z[0]-z[1]),
            0.25*(1.0+z[0])*(1.0+z[1])*(-1.0+z[0]+z[1]),
            0.25*(1.0-z[0])*(1.0+z[1])*(-1.0-z[0]+z[1]),
            0.5*(1.0-z[0]**2.0)*(1.0-z[1]),
            0.5*(1.0+z[0])*(1.0-z[1]**2.0),
            0.5*(1.0-z[0]**2.0)*(1.0+z[1]),
            0.5*(1.0-z[0])*(1.0-z[1]**2.0)
        ]).T

    def dpsis(self, z: np.ndarray) -> np.ndarray:
        """Calculates the shape functions derivatives of a given natural coordinates

        Args:
            z (np.ndarray): Natural coordinates matrix

        Returns:
            np.ndarray: Shape function derivatives evaluated in Z points
        """

        return np.array(
            [[-0.25*(z[1]-1.0)*(2.0*z[0]+z[1]), -0.25*(z[0]-1.0)*(2.0*z[1]+z[0])],
             [-0.25*(z[1]-1.0)*(2.0*z[0]-z[1]),
              0.25*(z[0]+1.0)*(2.0*z[1]-z[0])],
             [0.25*(z[1]+1.0)*(2.0*z[0]+z[1]),
              0.25*(z[0]+1.0)*(2.0*z[1]+z[0])],
             [0.25*(z[1]+1.0)*(2.0*z[0]-z[1]), -
              0.25*(z[0]-1.0)*(2.0*z[1]-z[0])],
             [(z[1]-1.0)*z[0], 0.5*(z[0]**2.0-1.0)],
             [-0.5*(z[1]**2.0-1.0), -z[1]*(z[0]+1.0)],
             [-(z[1]+1.0)*z[0], -0.5*(z[0]**2.0-1.0)],
             [0.5*(z[1]**2.0-1.0), z[1]*(z[0]-1.0)]])
