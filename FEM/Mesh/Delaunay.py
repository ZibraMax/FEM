"""Defines 2D domains by Delaunay triangulations
"""


import triangle as tr
import numpy as np
import copy
import matplotlib.pyplot as plt
from .Geometry import Geometry
from ..Utils import roundCorner, giveCoordsCircle


class Delaunay(Geometry):
    """Generate Delaunay triangulation using Triangle

    Args:
            vertices (list): matrix containing the domain vertices coordinates
            params (str): Triangulation parameters, use the aux function _strdelaunay
            nvn (int, optional): Number of variables per node. Defaults to 1.
            holes (list, optional): A list of holes. Defaults to None.
            fillets (list, optional): A list of fillets. Defaults to None.
    """

    def __init__(self, vertices: list, params: str, nvn: int = 1, holes_dict=None, fillets=None) -> None:
        """Generate Delaunay triangulation

        Args:
                vertices (list): matrix containing the domain vertices coordinates
                params (str): Triangulation parameters, use the aux function _strdelaunay
                nvn (int, optional): Number of variables per node. Defaults to 1.
                holes (list, optional): A list of holes. Defaults to None.
                fillets (list, optional): A list of fillets. Defaults to None.
        """
        mask = copy.deepcopy(vertices)
        # try:
        # mask = mask.tolist()
        # except:
        # pass
        seg = []
        for i in range(len(vertices)-1):
            seg.append([i, i+1])
        seg.append([i+1, 0])
        hh = []
        mascarita = copy.deepcopy(seg)
        if fillets:
            for fillet in fillets:
                S1 = seg[fillet['start_segment']]
                S2 = seg[fillet['end_segment']]
                for i, maskarita in enumerate(mascarita):
                    if maskarita == S1:
                        indice_importante = i

                mizq = mascarita[:indice_importante]
                mder = mascarita[indice_importante:]

                P1 = vertices[S1[0]]
                P2 = vertices[S2[1]]
                P = vertices[S1[1]]
                r = fillet['r']
                n = fillet['n']
                if not S1[1] == S2[0]:
                    raise Exception('The fillet segments are not valid')
                O, sa, a = roundCorner(P1, P2, P, r)
                f_vertices, f_segments = giveCoordsCircle(O, r, sa, a, n, True)
                vertices[S1[1]] = np.array(f_vertices[0]).tolist()
                sp = (np.array(f_segments)+len(vertices)-2)[1:].tolist()
                seg += [[S1[1], sp[1][0]]]+sp[1:]
                mder = mder[1:]
                spp = copy.deepcopy(sp)
                ss1 = copy.deepcopy(S1)
                if mder:
                    mder[0][0] = spp[-1][-1]
                mascarita = mizq+[[mizq[-1][-1], ss1[1]],
                                  [ss1[1], spp[1][0]]]+spp[1:]+mder
                vertices += np.array(f_vertices)[1:-1].tolist()
                seg[fillet['end_segment']][0] = len(vertices)-1
                # vertices += [O]

        original = dict(vertices=np.array(vertices), segments=np.array(seg))
        if holes_dict:
            for hole in holes_dict:
                hh += [hole['center']]
                seg += (np.array(hole['segments'])+len(vertices)).tolist()
                hole['segments'] = (
                    np.array(hole['segments'])+len(vertices)).tolist()
                vertices += np.array(hole['vertices']).tolist()
            original = dict(vertices=np.array(vertices),
                            segments=np.array(seg), holes=hh)
        triangular = tr.triangulate(original, params)
        dictionary = triangular['triangles'].tolist()
        tipos = np.zeros([len(dictionary)]).astype(str)
        if 'o2' in params:
            tipos[:] = 'T2V'
        else:
            tipos[:] = 'T1V'
        gdls = triangular['vertices'].tolist()
        if tipos[0] == 'T2V':
            for dicc in dictionary:
                a1 = dicc[5]
                a2 = dicc[3]
                a3 = dicc[4]
                dicc[3] = a1
                dicc[4] = a2
                dicc[5] = a3
        Geometry.__init__(self, dictionary, gdls, tipos,
                          nvn=nvn, segments=seg)
        mask = []
        for segmento in mascarita:
            mask += [gdls[segmento[0]]]
        if fillets:
            self.mask = mask

        self.holes = holes_dict
        self.fillets = fillets

    @ staticmethod
    def _strdelaunay(constrained: bool = True, delaunay: bool = True, a: float = None, q: float = None, o: int = 1) -> str:
        """Create a string for the delaunay triangulation constructor

        Args:
                constrained (bool, optional): Makes the triangulation constrained. Defaults to True.
                delaunay (bool, optional): Makes all triangles delaunay. Defaults to True.
                a (float, optional): Maximum area of triange. Defaults to None.
                q (float, optional): Minimum triangle angle <=35. Defaults to None.
                o (int, optional): Order of element if 2, quadratic elements are generated. Defaults to 1.

        Returns:
                str: A string containing the input parameters for the Delaunay1V constructor
        """
        p = ''
        if o == 2:
            o = '-o2'
        else:
            o = ''
        if constrained:
            p = 'p'
        if a == None:
            a = ''
        else:
            a = 'a'+format(a)
        D = ''
        if delaunay:
            D = 'D'
        if q == None:
            q = ''
        else:
            if type(q) == int:
                if q > 35:
                    raise "No se puede crear una triangulacion con angulos menores a 35 grados"
            q = 'q'+format(q)
        return p+a+D+q+'i'+o
