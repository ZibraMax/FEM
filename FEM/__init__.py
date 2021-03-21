"""FEM implementation
"""

from .__version__ import __version__
from .Elements import *
from .Core import *
from .Torsion2D import *
from .EDO1D import *
from .EulerBernoulliBeam import *
from .PlaneStress import *
from .PlaneStressNonLocal import *
from .PlaneStrain import *

__author__ = "Arturo Rodriguez - da.rodriguezh@uniandes.edu.co"

__all__ = ["Elements", "Core", "Torsion2D", "EDO1D", "EulerBernoulliBeam",
           "PlaneStress", "PlaneStrain", "PlaneStressNonLocal", "Mesh"]
