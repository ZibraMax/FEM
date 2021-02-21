import numpy as np

class LinearScheme():
	def __init__(this, n):
		this.Z,this.W = np.polynomial.legendre.leggauss(n)
		this.domain = np.array([-1] + this.Z.tolist() + [1])
	def integrate(this,f):
		integral = 0
		for w,z in zip(this.W,this.Z):
			integral += f(z)*w
		return integral