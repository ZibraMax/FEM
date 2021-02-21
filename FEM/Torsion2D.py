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
			dpx = []
			for i in range(len(_j)):
				dpx.append(_j[i] @ dpz.T[i].T)
			# dpx = (_j @ dpz.T).T #Shape function derivatives in global coordinates
			dpx = np.array(dpx).T
			for k in range(len(_x)): #Iterate over gauss points on domain
				for i in range(e.n): #This part must be vectorized
					for j in range(e.n):
						e.Ke[i,j] += (dpx[k][0][i]*dpx[k][0][j] + dpx[k][1][i]*dpx[k][1][j])*detjac[k]
					e.Fe[i][0] += 2*this.G*this._phi*_p[k][i]*detjac[k]

	def postProcess(this):
		X = []
		Y = []
		U = []
		fig = plt.figure()
		ax = fig.add_subplot(projection='3d')
		for e in this.elements:
			_x,_u=e.giveSolution()
			X+=_x[0].tolist()
			Y+=_x[1].tolist()
			U+=_u[0].tolist()
		surf = ax.plot_trisurf(X,Y,U,cmap='magma')
		surf._facecolors2d=surf._facecolors3d
		surf._edgecolors2d=surf._edgecolors3d
		# np.savetxt('U.csv',this.U,fmt='%s')
		# print(this.U)

