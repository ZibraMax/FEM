from ..Element import *
from .LinearScheme import *
import numpy as np
import matplotlib.pyplot as plt

class Element1D(Element,LinearScheme):
	def __init__(this,coords,gdl,n):
		coords = np.array(coords).reshape([len(coords),1])
		_coords = np.array([coords[0][0],coords[-1][0]])
		Element.__init__(this,coords,_coords,gdl)
		LinearScheme.__init__(this,n)

	def draw(this,m=100):
		_z = np.linspace(*this.domain.reshape([2]),m)
		_x,_p = this.T(_z)
		_y = 0
		for i in range(this.n):
			plt.plot(_x,_p[i],'--',label=r'$\psi_{'+format(i)+r'}$')
			_y+=_p[i]
		plt.plot(_x,[0]*len(_x),'-',color='black',label='Element')
		plt.plot(this.coords.T[0],[0]*len(this.coords.T[0]),'o',color='blue',label='Nodes')
		plt.plot(_x,_y,'-.',label=r'$\sum_{i=0}^{n}\psi$')
		plt.legend()
		plt.grid()
		
