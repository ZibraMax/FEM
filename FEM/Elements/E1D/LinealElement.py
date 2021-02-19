from .Element1D import *
class LinealElement(Element1D):
	def __init__(this,coords,gdl,n=2):
		Element1D.__init__(this,coords,gdl,n)
	def psis(this,z):
		return np.array([0.5*(1.0-z),0.5*(1.0+z)])
	def dpsis(this,z):
		return np.array([[-0.5,0.5]])