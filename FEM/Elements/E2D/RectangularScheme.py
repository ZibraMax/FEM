import numpy as np

class RectangularScheme():
	def __init__(this, n):
		DENSIDAD = 10
		XX,YY = np.meshgrid(np.linspace(-1,1,DENSIDAD), np.linspace(-1,1,DENSIDAD))
		_z = XX.reshape([DENSIDAD**2,1])[:,0]
		_n = YY.reshape([DENSIDAD**2,1])[:,0]
		this.domain = np.array([_z,_n]).T
		_Z,_W = np.polynomial.legendre.leggauss(n)
		this.Z = []
		this.W = []
		for i,z in enumerate(_Z):
			for j,n in enumerate(_Z):
				this.Z.append([z,n])
				this.W.append(_W[i]*_W[j])

		this.Z = np.array(this.Z)
		this.W = np.array(this.W)