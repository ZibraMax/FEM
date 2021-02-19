import numpy as np

class LinearScheme():
	def __init__(this, n):
		this.domain = np.array([[-1.0],[1.0]])
		this.Z,this.W = np.polynomial.legendre.leggauss(n)
	def integrate(this,f):
		integral = 0
		for w,z in zip(this.W,this.Z):
			integral += f(z)*w
		return integral