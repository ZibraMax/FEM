from tqdm import tqdm
import numpy as np
class Core():
	def __init__(this,geometry):
		this.geometry = geometry
		this.ngdl = this.geometry.ngdl
		this.K = np.zeros([this.ngdl,this.ngdl])
		this.F = np.zeros([this.ngdl,1])
		this.Q = np.zeros([this.ngdl,1])
		this.U = np.zeros([this.ngdl,1])
		this.S = np.zeros([this.ngdl,1])
		this._K = np.zeros([this.ngdl,this.ngdl])
		this.cbe = this.geometry.cbe
		this.cbn = this.geometry.cbn
		this.elements = this.geometry.elements

	def ensembling(this):
		print('Ensembling equation system...')
		for e in tqdm(this.elements,unit='Element'):
			this.K[np.ix_(e.gdl.flatten(),e.gdl.flatten())] += e.Ke
			this.F[np.ix_(e.gdl.flatten())] += e.Fe
			this.Q[np.ix_(e.gdl.flatten())] += e.Qe
		this._K = np.copy(this.K)
		print('Done!')
	def borderConditions(this):
		for i in this.cbn:
			this.Q[int(i[0])] = i[1]
		for i in this.cbe:
			ui = np.zeros([this.ngdl, 1])
			ui[int(i[0])] = i[1]
			vv = np.dot(this.K, ui)
			this.S = this.S - vv
			this.K[int(i[0]), :] = 0
			this.K[:, int(i[0])] = 0
			this.K[int(i[0]), int(i[0])] = 1
		this.S = this.S + this.F + this.Q
		for i in this.cbe:
			this.S[int(i[0])] = i[1]
	def solveES(this,path=''):
		print('Solving equation system...')
		this.U = np.linalg.solve(this.K,this.S)
		if not path == '':
			np.savetxt(path,this.U,delimiter=',')
		for e in this.elements:
			e.setUe(this.U)
		print('Done!')
	def solve(this):
		this.elementMatrices()
		this.ensembling()
		this.borderConditions()
		this.solveES()
		this.postProcess()

	def elementMatrices(this): #This methods must be implemented in children classes
		pass
	def postProcess(this): #This methods must be implemented in children classes
		pass