import numpy as np
class Element():
	def __init__(this,coords,_coords,gdl):
		this.coords = coords
		this._coords = _coords
		this.gdl = gdl
		this.n = len(this.gdl)

	def T(this,z):
		x = 0
		p = this.psis(z)
		if len(this.coords[0])==1:
			for i in range(len(p)):
				x += this.coords[i].T * p[i]
		else:
			x = [0]*len(this.coords[0])
			for i in range(len(p)):
				for j in range(len(this.coords[i])):
					x[j] += this.coords[i][j] * p[i]
			x = np.array(x)
		return x,p

	def J(this, z):
		return this.coords.T @ this.dpsis(z).T
