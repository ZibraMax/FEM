from .Element2D import *
from .TriangularScheme import *

class LTriangular(Element2D,TriangularScheme):
	def __init__(this,coords,gdl,n=2):
		coords = np.array(coords)
		Element2D.__init__(this,coords,coords,gdl)
		TriangularScheme.__init__(this,n)
	def psis(this,z):
		return np.array([
			1.0-z[0]-z[1],
			z[0],
			z[1]])
	def dpsis(this,z):
		kernell = (z[0]-z[0])
		return np.array(
			[[-1.0*(1+kernell),1.0*(1+kernell),0.0*(1+kernell)],
			[-1.0*(1+kernell),0.0*(1+kernell),1.0*(1+kernell)]])
