from .Element2D import *
from .RectangularScheme import *

class Quadrilateral(Element2D,RectangularScheme):
	def __init__(this,coords,gdl,n=2):
		coords = np.array(coords)
		Element2D.__init__(this,coords,coords,gdl)
		RectangularScheme.__init__(this,n)
	def psis(this,z):
		return np.array(
			[0.25*(1.0-z[0])*(1.0-z[1]),
			0.25*(1.0+z[0])*(1.0-z[1]),
			0.25*(1.0+z[0])*(1.0+z[1]),
			0.25*(1.0-z[0])*(1.0+z[1])])
	def dpsis(this,z):
		return np.array(
			[[0.25*(z[1]-1.0),
			-0.25*(z[1]-1.0),
			0.25*(z[1]+1.0),
			-0.25*(1.0+z[1])],
			[0.25*(z[0]-1.0),
			-0.25*(z[0]+1.0),
			0.25*(1.0+z[0]),
			0.25*(1.0-z[0])]])

