from tqdm import tqdm
import numpy as np
class Core():
	def __init__(self,geometry):
		self.geometry = geometry
		self.ngdl = self.geometry.ngdl
		self.K = np.zeros([self.ngdl,self.ngdl])
		self.F = np.zeros([self.ngdl,1])
		self.Q = np.zeros([self.ngdl,1])
		self.U = np.zeros([self.ngdl,1])
		self.S = np.zeros([self.ngdl,1])
		self._K = np.zeros([self.ngdl,self.ngdl])
		self.cbe = self.geometry.cbe
		self.cbn = self.geometry.cbn
		self.elements = self.geometry.elements

	def ensembling(self):
		print('Ensembling equation system...')
		for e in tqdm(self.elements,unit='Element'):
			self.K[np.ix_(e.gdlm,e.gdlm)] += e.Ke
			self.F[np.ix_(e.gdlm)] += e.Fe
			self.Q[np.ix_(e.gdlm)] += e.Qe
		self._K = np.copy(self.K)
		print('Done!')
	def borderConditions(self):
		print('Border conditions...')
		for i in tqdm(self.cbn,unit=' Natural'):
			self.Q[int(i[0])] = i[1]
		for i in tqdm(self.cbe,unit=' Essential'):
			ui = np.zeros([self.ngdl, 1])
			ui[int(i[0])] = i[1]
			vv = np.dot(self.K, ui)
			self.S = self.S - vv
			self.K[int(i[0]), :] = 0
			self.K[:, int(i[0])] = 0
			self.K[int(i[0]), int(i[0])] = 1
		self.S = self.S + self.F + self.Q
		for i in self.cbe:
			self.S[int(i[0])] = i[1]
	def solveES(self,path=''):
		print('Solving equation system...')
		self.U = np.linalg.solve(self.K,self.S)
		if not path == '':
			np.savetxt(path,self.U,delimiter=',')
		for e in self.elements:
			e.setUe(self.U)
		print('Done!')
	def solve(self,path=''):
		self.elementMatrices()
		self.ensembling()
		self.borderConditions()
		self.solveES(path)
		self.postProcess()
	
	def solveFromFile(self,file):
		# self.elementMatrices()
		# self.ensembling()
		# self.borderConditions()
		self.U = np.loadtxt(file)
		for e in self.elements:
			e.setUe(self.U)
		self.postProcess()

	def elementMatrices(self): #self methods must be implemented in children classes
		pass
	def postProcess(self): #self methods must be implemented in children classes
		pass