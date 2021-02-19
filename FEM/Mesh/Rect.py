import numpy as np
from .Geometria import *

class Rect(Geometria):
	def __init__(this,file):
		gdls = []
		diccionarios = []
		diccionariosnl = []
		vertices = [-1]
		segmentos = []
		with open(file) as archivo:
			params = archivo.readline().split('\t')
			ngdls = int(params[0])
			nele = int(params[1])
			tipos = np.zeros([nele]).astype(str)
			tipos[:] = 'C1V'
			for _ in range(ngdls):
				linea = archivo.readline().split('\t')
				gdls.append([float(linea[0]),float(linea[1])])
			for i in range(nele):
				linea = list(map(lambda x: int(x)-1,archivo.readline().split('\t')))
				diccionarios.append(linea)
				if len(linea)==3:
					tipos[i]='T1V'
				if len(linea)==6:
					tipos[i]='T2V'
				if len(linea)==4:
					tipos[i]='C1V'
				if len(linea)==8:
					tipos[i]='C2V'
			for _ in range(nele):
				linea = list(map(lambda x: int(x)-1,archivo.readline().split('\t')[1:]))
				diccionariosnl.append(linea)
		super().__init__(vertices, diccionarios, gdls, tipos)
		this.diccionariosnl = diccionariosnl
		this.cbe = []

	def generarCBdesdeBordeX(this, borde, valor=0):
		cb = []
		nodos = this.darNodosCB(borde)
		cbe = np.zeros([len(nodos), 2])
		cbe[:, 0] = nodos*2
		cbe[:, 1] = valor
		cb += cbe.tolist()
		return cb

	def generarCBdesdeBordeY(this, borde, valor=0):
		cb = []
		nodos = this.darNodosCB(borde)
		cbe = np.zeros([len(nodos), 2])
		cbe[:, 0] = nodos*2+1
		cbe[:, 1] = valor
		cb += cbe.tolist()
		return cb

	def generarCBXdesdeCoordenada(this,x,y,valor=0):
		masCercano1 = None
		d1 = np.Inf
		for i,gdl in enumerate(this.gdls):
			r1 = np.sqrt((x-gdl[0])**2+(y-gdl[1])**2)
			if r1 < d1:
				d1 = r1
				masCercano1 = i
		return [[i*2,valor]]

	def generarCBYdesdeCoordenada(this,x,y,valor=0):
		masCercano1 = None
		d1 = np.Inf
		for i,gdl in enumerate(this.gdls):
			r1 = np.sqrt((x-gdl[0])**2+(y-gdl[1])**2)
			if r1 < d1:
				d1 = r1
				masCercano1 = i
		return [[i*2+1,valor]]

	def generarCBdesdeBorde(this, borde, valor=[0,0]):
		return this.generarCBdesdeBordeX(borde, valor[0])+this.generarCBdesdeBordeY(borde, valor[1])