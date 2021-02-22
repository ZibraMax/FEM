import numpy as np
class Element():
	def __init__(this,coords,_coords,gdl):
		this.coords = coords
		this._coords = _coords
		this.gdl = gdl
		this.n = int(len(this.gdl)*len(this.gdl[0]))
		this.Ke = np.zeros([this.n,this.n])
		this.Fe = np.zeros([this.n,1])
		this.Ue = np.zeros(this.gdl.shape)
		this.Qe = np.zeros([this.n,1])

	def T(this,z):
		p = this.psis(z)
		return p@this.coords,p

	def inverseMapping(this,x,n=100):
		x = np.array(x)
		tol = 1*10**(-6)
		zeta = np.zeros(x.shape)+0.1
		for i in range(n):
			xi = x - this.T(zeta)[0]
			J = np.linalg.inv(this.J(zeta)[0])
			so = list(xi.shape)
			xi = xi.reshape(list(xi.shape)+[1])
			dz = J@xi
			zeta += dz.reshape(so)
			if np.max(np.abs(dz))<tol:
				break
		return zeta
	def J(this, z):
		dpsis = this.dpsis(z).T
		return dpsis @ this.coords, dpsis

	def giveSolution(this,SVSolution=False):
		_z = this.domain
		_x,_p = this.T(_z.T)
		if SVSolution:
			j,dpz = this.J(_z.T) #TODO Revisar con Reddy
			dpx =  np.linalg.inv(j) @ dpz
			# print((this.Ue @ np.transpose(dpx,axes=[0,2,1])).shape)
			return _x,this.Ue@_p.T, this.Ue @ np.transpose(dpx,axes=[0,2,1]) #TODO REVISAR VS
		return _x,this.Ue@_p.T

	def setUe(this,U):
		for i in range(len(this.gdl)):
			this.Ue[i] = U[np.ix_(this.gdl[i])].flatten()

	def integrate(this,f):
		integral = 0
		for w,z in zip(this.W,this.Z):
			integral += f(z)*w
		return integral

