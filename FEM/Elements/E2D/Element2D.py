from ..Element import *
import numpy as np
import matplotlib.pyplot as plt
import numpy as np

class Element2D(Element):
	def __init__(this,coords,_coords,gdl):
		Element.__init__(this,coords,_coords,gdl)

	def draw(this):
		_z = this.domain
		_x,_p = this.T(_z.T)
		fig = plt.figure()
		ax = fig.add_subplot(projection='3d')
		l = []
		l.append('Element')
		l.append('Nodes')
		for i in range(this.n):
			surf = ax.plot_trisurf(*_x,_p[i],color=(0,0,0,0), edgecolor=(np.random.random(),np.random.random(),np.random.random()))
			surf._facecolors2d=surf._facecolors3d
			surf._edgecolors2d=surf._edgecolors3d
			l.append(r'$\psi_{'+format(i)+r'}$')
		__coords = np.array(this._coords.tolist()+[this._coords[0].tolist()]).T
		ax.plot(*__coords,[0]*len(__coords.T),'-',color='black')
		ax.plot(*this.coords.T,[0]*len(this.coords),'o',color='blue')
		ax.legend(l)
	def jacobianGraph(this):
		_z = this.domain
		_x,_p = this.T(_z.T)
		_j = this.J(_z.T)
		__j = np.linalg.det(_j)
		fig = plt.figure()
		ax = fig.add_subplot(projection='3d')
		l = []
		surf = ax.plot_trisurf(*_x,__j,cmap='magma')
		surf._facecolors2d=surf._facecolors3d
		surf._edgecolors2d=surf._edgecolors3d
		l.append('Element')
		l.append('Nodes')
		l.append(r'$|J|$')
		cbar = fig.colorbar(surf)
		__coords = np.array(this._coords.tolist()+[this._coords[0].tolist()]).T
		ax.plot(*__coords,[0]*len(__coords.T),'-',color='black')
		ax.plot(*this.coords.T,[0]*len(this.coords),'o',color='blue')
		ax.legend(l)


