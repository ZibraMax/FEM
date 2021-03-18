import triangle as tr

import numpy as np
import matplotlib.pyplot as plt

from .Geometry import *


class Delaunay1V(Geometry):
    def __init__(self, vertices, params, nvn=1):
        """Generate Delaunay triangulation

        Args:
                vertices (list): matrix containing the domain vertices coordinates
                params (str): Triangulation parameters, use the aux function _strdelaunay
                nvn (int, optional): Number of variables per node. Defaults to 1.
        """
        seg = []
        for i in range(len(vertices)-1):
            seg.append([i, i+1])
        seg.append([i+1, 0])
        original = dict(vertices=np.array(vertices), segments=np.array(seg))
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
        Geometry.__init__(self, dictionary, gdls, tipos, nvn=nvn, segments=seg)
        self.mask = vertices


def _strdelaunay(constrained=True, delaunay=True, a=None, q=None, o=1):
    """Create a string for the delaunay triangulation constructor

    Args:
            constrained (bool, optional): Makes the triangulation constrained. Defaults to True.
            delaunay (bool, optional): Makes all triangles delaunay. Defaults to True.
            a ([type], optional): Maximum area of triange. Defaults to None.
            q ([type], optional): Minimum triangle angle <=35. Defaults to None.
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
