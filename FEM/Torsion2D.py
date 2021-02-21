from .Core import *
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

class Torsion2D(Core):
	def __init__(this,geometry,G,phi):
		this.G = G
		this._phi = phi
		geometry.cbeAllBorders(0)
		Core.__init__(this,geometry)

	def elementMatrices(this):
		for e in tqdm(this.elements,unit='Element'):
			_x,_p = e.T(e.Z.T) #Gauss points in global coordinates and Shape functions evaluated in gauss points
			j,dpz = e.J(e.Z.T) #Jacobian evaluated in gauss points and shape functions derivatives in natural coordinates
			detjac = np.linalg.det(j)
			_j = np.linalg.inv(j) #Jacobian inverse
			dpx = dpz @ _j #Shape function derivatives in global coordinates
			for k in range(len(_x)): #Iterate over gauss points on domain
				for i in range(e.n): #This part must be vectorized
					for j in range(e.n):
						e.Ke[i,j] += (dpx[k][i][0]*dpx[k][j][0] + dpx[k][i][1]*dpx[k][j][1])*detjac[k]
					e.Fe[i][0] += 2*this.G*this._phi*_p[k][i]*detjac[k]

	def postProcess(this):
		X = []
		Y = []
		U1 = []
		U2 = []
		U3 = []
		U4 = []
		fig = plt.figure()
		ax1 = fig.add_subplot(2,2,1)
		ax2 = fig.add_subplot(2,2,2)
		ax3 = fig.add_subplot(2,2,3)
		ax4 = fig.add_subplot(2,2,4)
		for e in this.elements:
			_x,_u,du=e.giveSolution(True)
			X+=_x.T[0].tolist()
			Y+=_x.T[1].tolist()
			U2+=du[0,0].tolist()
			U3+=du[1,0].tolist()
			U1+=_u[0].tolist()
			U4+=np.sqrt(du[0,0]**2 + du[1,0]**2).tolist()
		surf = ax1.tricontourf(X,Y,U1,cmap='magma')
		cbar = fig.colorbar(surf,ax=ax1)
		surf = ax2.tricontourf(X,Y,U2,cmap='magma')
		cbar = fig.colorbar(surf,ax=ax2)
		surf = ax3.tricontourf(X,Y,U3,cmap='magma')
		cbar = fig.colorbar(surf,ax=ax3)
		surf = ax4.tricontourf(X,Y,U4,cmap='magma')
		cbar = fig.colorbar(surf,ax=ax4)

