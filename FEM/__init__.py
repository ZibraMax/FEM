"""FEM implementation for N dimensional and M variables per node.

"""
__author__ = "Arturo Rodriguez - da.rodriguezh@uniandes.edu.co"

from .__version__ import __version__
from .Elements import *
from .Mesh import *
from .Core import *
from .Torsion2D import *
from .EDO1D import *
from .EulerBernoulliBeam import *
from .PlaneStress import *
from .PlaneStressNonLocal import *
from .PlaneStrain import *
from .Heat1D import *
from .Heat2D import *
