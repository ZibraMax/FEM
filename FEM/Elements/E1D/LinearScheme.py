import numpy as np

class LinearScheme():
	def __init__(self, n):
		self.Z,self.W = np.polynomial.legendre.leggauss(n)
		self.domain = self.Z
	def integrate(self,f):
		integral = 0
		for w,z in zip(self.W,self.Z):
			integral += f(z)*w
		return integral