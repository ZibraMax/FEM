from .Element1D import *
class QuadraticElement(Element1D):
	def __init__(this,coords,gdl,n=2):
		Element1D.__init__(this,coords,gdl,n)
	def psis(this,z):
		zm1 = z+1.0
		return np.array([1.0-3.0/2.0*zm1+zm1*zm1/2.0,2.0*zm1*(1.0-zm1/2.0),z/2.0*zm1]).T
	def dpsis(this,z):
		return np.array([[z-0.5],[-2.0*z],[z+0.5]])