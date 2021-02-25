from .Core import *
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

class PlaneStress(Core):
	
	def __init__(self,geometry,E,v,t):
		
		if type(t)==float:
			t = [t]*len(geometry.elements)
		if type(E)==float:
			E = [E]*len(geometry.elements)
		if type(v)==float:
			v = [v]*len(geometry.elements)
		self.t = t
		self.E = E
		self.v = v
		self.C11 = []
		self.C12 = []
		self.C66 = []

		for i in range(len(self.E)):
			C11 = self.E[i] / (1 - self.v[i] ** 2)
			C12 = self.v[i] * self.E[i] / (1 - self.v[i] ** 2)
			C66 = self.E[i] / 2 / (1 + self.v[i])
			self.C11.append(C11)
			self.C12.append(C12)
			self.C66.append(C66)
		if geometry.nvn == 1:
			print('Border conditions lost, please usea a geometry with 2 variables per node (nvn=2)')
			geometry.nvn = 2
			geometry.cbe = []
			geometry.cbn = []
			geometry.initialize()
		Core.__init__(self,geometry)

	def elementMatrices(self):

# EKuu = lambda z, n: (C11 * dfdx(z, n, i) * dfdx(z, n, j) + C66 * dfdy(z, n, i) * dfdy(z, n, j)) * np.linalg.det(J(z, n))
# EKuv = lambda z, n: (C12 * dfdx(z, n, i) * dfdy(z, n, j) + C66 * dfdy(z, n, i) * dfdx(z, n, j)) * np.linalg.det(J(z, n))
# EKvu = lambda z, n: (C12 * dfdy(z, n, i) * dfdx(z, n, j) + C66 * dfdx(z, n, i) * dfdy(z, n, j)) * np.linalg.det(J(z, n))
# EKvv = lambda z, n: (C11 * dfdy(z, n, i) * dfdy(z, n, j) + C66 * dfdx(z, n, i) * dfdx(z, n,j)) * np.linalg.det(J(z, n))
		ee = 0
		for e in tqdm(self.elements,unit='Element'):
			m = len(e.gdl.T)
			Kuu = np.zeros([m,m])
			Kuv = np.zeros([m,m])
			Kvu = np.zeros([m,m])
			Kvv = np.zeros([m,m])
			Fu = np.zeros([m,1])
			Fv = np.zeros([m,1])
			_x,_p = e.T(e.Z.T) #Gauss points in global coordinates and Shape functions evaluated in gauss points
			jac,dpz = e.J(e.Z.T) #Jacobian evaluated in gauss points and shape functions derivatives in natural coordinates
			detjac = np.linalg.det(jac)
			_j = np.linalg.inv(jac) #Jacobian inverse
			dpx = _j @ dpz #Shape function derivatives in global coordinates
			for i in range(m): #self part must be vectorized
				for j in range(m):
					for k in range(len(e.Z)): #Iterate over gauss points on domain
						Kuu[i,j] += self.t[ee]*(self.C11[ee]*dpx[k,0,i]*dpx[k,0,j]+self.C66[ee]*dpx[k,1,i]*dpx[k,1,j])*detjac[k]*e.W[k]
						Kuv[i,j] += self.t[ee]*(self.C12[ee]*dpx[k,0,i]*dpx[k,1,j]+self.C66[ee]*dpx[k,1,i]*dpx[k,0,j])*detjac[k]*e.W[k]
						Kvu[i,j] += self.t[ee]*(self.C12[ee]*dpx[k,1,i]*dpx[k,0,j]+self.C66[ee]*dpx[k,0,i]*dpx[k,1,j])*detjac[k]*e.W[k]
						Kvv[i,j] += self.t[ee]*(self.C11[ee]*dpx[k,1,i]*dpx[k,1,j]+self.C66[ee]*dpx[k,0,i]*dpx[k,0,j])*detjac[k]*e.W[k]
				for k in range(len(e.Z)): #Iterate over gauss points on domain
					e.Fe[i][0] += 0
			subm = np.linspace(0, 2*m-1,2*m).reshape([2, m]).astype(int)
			e.Ke[np.ix_(subm[0],subm[0])] += Kuu
			e.Ke[np.ix_(subm[0],subm[1])] += Kuv
			e.Ke[np.ix_(subm[1],subm[0])] += Kvu
			e.Ke[np.ix_(subm[1],subm[1])] += Kvv
			ee+=1
			# e.Fe[:,0] = 2*self.G*self._phi*detjac@_p
			# e.Ke = (np.transpose(dpx,axes=[0,2,1]) @ dpx).T @ detjac
	def postProcess(self,mult=800):
		print('Post-processing...')
		X = []
		Y = []
		U1 = []
		U2 = []
		U3 = []
		U4 = []
		Ux = []
		Uy = []
		Ux0 = []
		Uy0 = []
		fig = plt.figure()

		gs = gridspec.GridSpec(4, 2)

		ax1 = fig.add_subplot(gs[0,0])
		ax2 = fig.add_subplot(gs[0,1])
		ax3 = fig.add_subplot(gs[1,0])
		ax4 = fig.add_subplot(gs[1,1])
		ax5 = fig.add_subplot(gs[2:,:])
		for e in tqdm(self.elements,unit='Element'):
			_x,_u,du=e.giveSolution(True)
			X+=_x.T[0].tolist()
			Y+=_x.T[1].tolist()
			U1+=(du[:,0,0]).tolist()
			U2+=(du[:,0,1]).tolist()
			U3+=(du[:,1,0]).tolist()
			U4+=(du[:,0,1]).tolist()
			coordsNuevas = e.coords + e.Ue.T * mult
			Ux += coordsNuevas.T.tolist()[0]
			Uy += coordsNuevas.T.tolist()[1]
			Ux0 += e.coords.T.tolist()[0]
			Uy0 += e.coords.T.tolist()[1]

		surf = ax1.tricontourf(X,Y,U1,cmap='magma')
		plt.colorbar(surf,ax=ax1)
		ax1.set_title(r'$\frac{dU}{dx}$')

		surf = ax2.tricontourf(X,Y,U2,cmap='magma')
		plt.colorbar(surf,ax=ax2)
		ax2.set_title(r'$\frac{dU}{dy}$')

		surf = ax3.tricontourf(X,Y,U3,cmap='magma')
		plt.colorbar(surf,ax=ax3)
		ax3.set_title(r'$\frac{dV}{dx}$')

		surf = ax4.tricontourf(X,Y,U4,cmap='magma')
		plt.colorbar(surf,ax=ax4)
		ax4.set_title(r'$\frac{dV}{dy}$')

		ax5.plot(Ux0,Uy0,'--',color='gray',alpha=0.7)
		ax5.plot(Ux,Uy,'-',color='black')
		print('Done!')
