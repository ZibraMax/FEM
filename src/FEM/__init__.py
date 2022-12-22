"""FEM implementation for N dimensional and M variables per node.

"""
__author__ = "Arturo Rodriguez - da.rodriguezh@uniandes.edu.co"
__version__ = "1.0.34"
from .Elements import *
from .Geometry import *
from .Core import *
from .Torsion2D import *
from .EDO1D import *
from .EulerBernoulliBeam import *
from .Elasticity2D import *
from .NonLinealExample import *
from .Heat1D import *
from .Heat2D import *

from .Elasticity3D import *

import sys
import inspect


def importJSON(json_file, aditional, **kargs):
    logging = FEMLogger()
    y = json.load(open(json_file))
    problemClass = y['properties']['problem']
    current_module = sys.modules[__name__]
    baseClass = current_module.__dict__[problemClass]
    coords = y['nodes']
    dcc = y['dictionary']
    types = y['types']
    regions = y['regions']
    ebc = y['ebc']
    nbc = y['nbc']
    nvn = y['nvn']
    solutions = y['solutions']
    sol = []
    sol_info = []

    for solution in solutions:
        sol += [np.array(solution['U'])]
        sol_info += [solution['info']]
    properties = y['properties']
    ndim = len(coords[0])
    if ndim == 1:
        geo = Geometry1D(dcc, coords, types, nvn, **kargs)
    elif ndim == 2:
        geo = Geometry2D(dcc, coords, types, nvn, **kargs)
    elif ndim == 3:
        geo = Geometry3D(dcc, coords, types, nvn, **kargs)
    else:
        raise Exception('Not valid geometry')
    params = inspect.signature(baseClass).parameters.keys()
    properties['geometry'] = geo
    for p in params:
        if not p in properties and not p == 'kargs':
            properties[p] = 0.5
            if p in aditional:
                properties[p] = aditional[p]
    for_Deletion = []
    for k in properties:
        if not k in params:
            for_Deletion += [k]
    for k in for_Deletion:
        properties.pop(k)

    O = baseClass(**properties)
    O.solver.solutions = sol
    O.solver.solutions_info = sol_info
    O.solver.setSolution(-1, elements=True)
    O.cbe = ebc
    O.cbn = nbc
    return O
