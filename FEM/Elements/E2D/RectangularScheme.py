import numpy as np

class RectangularScheme():
	def __init__(this, n):
		DENSIDAD = 10
		XX,YY = np.meshgrid(np.linspace(-1,1,DENSIDAD), np.linspace(-1,1,DENSIDAD))
		_z = XX.reshape([DENSIDAD**2,1])[:,0]
		_n = YY.reshape([DENSIDAD**2,1])[:,0]
		this.domain = np.array([_z,_n]).T
		this.Z,this.W = np.polynomial.legendre.leggauss(n)
		
	def integrate(this,f):
		integral = 0
		for wz,z in zip(this.W,this.Z):
			for wn,n in zip(this.W,this.Z):
				integral += f(np.array([z,n]))*wz*wn
		return integral