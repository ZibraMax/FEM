import numpy as np
from ..Utils import isBetween
import matplotlib.pyplot as plt
from ..Elements import *
from ..Elements.E1D import *
from .Geometry import *
import re

class Geometry1D(Geometry):
	def __init__(this, dictionary, gdls, types, nvn=1):
		this.nvn = nvn
		this.dictionary = dictionary
		this.elements = []
		this.gdls = gdls
		this.types = types
		this.cbe = []
		this.cbn = []
		this.ngdl = int(len(this.gdls)*this.nvn)
		this.generateElements()
	@staticmethod
	def loadmsh(filename):
		f = open(filename,'r')
		dicc = []
		gdls = []
		types = []
		seg = []
		cbe = []
		cbn = []
		nvn = 1
		p = list(map(int,f.readline().split('\t')))
		#[len(this.gdls),len(this.dictionary),len(this.segments),len(this.cbe),len(this.cbn),this.nvn]
		for _ in range(p[0]):
			gdls += [list(map(float,f.readline().split('\t')))]
		for _ in range(p[1]):
			types += [f.readline().split('\n')[0]]
		for _ in range(p[1]):
			dicc += [list(map(int,f.readline().split('\t')))]
		for _ in range(p[2]):
			seg += [list(map(int,f.readline().split('\t')))]
		for _ in range(p[3]):
			cbe += [list(map(int,f.readline().split('\t')))]
		for _ in range(p[4]):
			cbn += [list(map(int,f.readline().split('\t')))]
		nvn = p[5]
		f.close()
		print('File '+ filename + ' loaded')
		o = Geometry(dicc,gdls,types,nvn,seg)
		o.cbe = cbe
		o.cbn = cbn
		return o

	def generateElements(this):
		for i,d in enumerate(this.dictionary):
			coords = np.array(this.gdls)[np.ix_(d)]
			gdl = np.zeros([this.nvn,len(d)])
			for i in range(this.nvn):
				gdl[i,:] = (np.array(d)*this.nvn+i)
			gdl = gdl.astype(int)
			if this.types[i]=='L1V':
				element = LinealElement(coords,gdl)
			elif this.types[i]=='L2V':
				element = QuadraticElement(coords,gdl)
			this.elements.append(element)

	def show(this,texto=10,bolita=0,figsize=[17,10]):
		fig = plt.figure(figsize=figsize)
		ax = fig.add_subplot()

		ax.axes.set_aspect('equal')

		for i, e in enumerate(this.elements):
			coords = e._coords
			coords = np.array(coords.tolist() + [coords[0].tolist()])
			X = coords[:, 0]
			Y = coords[:, 1]
			ax.plot(X, Y, 'o-', color='black', zorder=-10)
			cx = this.centroids[i][0]
			cy = this.centroids[i][1]
			ax.plot(cx, cy, 'o', markersize=texto + bolita, color='yellow')
			ax.annotate(format(i), [cx, cy], size=texto, textcoords="offset points", xytext=(-0, -2.5), ha='center')
		try:
			verts = this.gdls
			segs = this.segments
			for i,seg in enumerate(segs):
				x0, y0 = verts[int(seg[0])]
				x1, y1 = verts[int(seg[1])]

				ax.fill(
					[x0, x1],
					[y0, y1],
					facecolor='none',
					edgecolor='b',
					linewidth=3,
					zorder=0,
				)
				cx = (x0+x1)*0.5
				cy = (y0+y1)*0.5
				ax.plot(cx, cy, 'o', markersize=texto + bolita, color='pink')
				ax.annotate(format(i), [cx, cy], size=texto, textcoords="offset points", xytext=(-0, -2.5), ha='center')
		except:
			pass
		ax.set_xlabel('x')
		ax.set_ylabel('y')
		ax.set_title('Domain')

		gdls = np.array(this.gdls)

		labels = np.linspace(0, gdls.shape[0] - 1, gdls.shape[0]).astype(int)

		ax.plot(gdls[:, 0], gdls[:, 1], 'o', markersize=texto+bolita, color='gray')

		for p, l in zip(gdls, labels):
			ax.annotate(l, p, size=texto, textcoords="offset points", xytext=(-0, -2.5), ha='center')

	def saveMesh(this,ProjectName):
		filename = ProjectName + '.msh'
		f = open(filename,'w')
		p = [len(this.gdls),len(this.dictionary),len(this.segments),len(this.cbe),len(this.cbn),this.nvn]
		f.write('\t'.join(list(map(str,p))) + '\n')
		for e in this.gdls:
			f.write('\t'.join(list(map(str,e)))+ '\n')
		for e in this.types:
			f.write(e+ '\n')
		for e in this.dictionary:
			f.write('\t'.join(list(map(str,e)))+ '\n')
		for e in this.segments:
			f.write('\t'.join(list(map(str,e)))+ '\n')
		for e in this.cbe:
			f.write('\t'.join(list(map(str,e)))+ '\n')
		for e in this.cbn:
			f.write('\t'.join(list(map(str,e)))+ '\n')
		f.close()
		print('File '+ filename + ' saved')

	def centroidsAndAreas(this):
		for i, e in enumerate(this.elements):
			coords = e._coords
			coords = np.array(coords.tolist() + [coords[0].tolist()])
			area = 0
			cx = 0
			cy = 0
			for j in range(len(coords)-1):
				area += coords[j][0]*coords[j+1][1]-coords[j+1][0]*coords[j][1]
				mult = (coords[j][0]*coords[j+1][1]-coords[j+1][0]*coords[j][1])
				cx += (coords[j][0]+coords[j+1][0])*mult
				cy += (coords[j][1]+coords[j+1][1])*mult
			this.areas.append(np.abs(area/2))
			this.centroids.append([cx/3/area,cy/3/area])

	def generateSegmentsFromCoords(this,p0,p1):
		masCercano1 = None
		d1 = np.Inf
		masCercano2 = None
		d2 = np.Inf
		for i,gdl in enumerate(this.gdls):
			r1 = np.sqrt((p0[0]-gdl[0])**2+(p0[1]-gdl[1])**2)
			r2 = np.sqrt((p1[0]-gdl[0])**2+(p1[1]-gdl[1])**2)
			if r1 < d1:
				d1 = r1
				masCercano1 = i
			if r2 < d2:
				d2 = r2
				masCercano2 = i
		this.segments.append([masCercano1,masCercano2])

	def generateBCFromCoords(this,x,y,value=0):
		masCercano1 = None
		d1 = np.Inf
		for i,gdl in enumerate(this.gdls):
			r1 = np.sqrt((x-gdl[0])**2+(y-gdl[1])**2)
			if r1 < d1:
				d1 = r1
				masCercano1 = i
		return [[i,valor]]

	def giveNodesOfSegment(this,segment,tol):
		a = []
		ps = np.array(this.gdls)[this.segments[segment]].tolist()
		for i, p in enumerate(this.gdls):
			if isBetween(ps[0], ps[1], p,tol):
				a.append(i)
		return np.array(a)

	def cbFromSegment(this,segment,value,nv=1,tol=1*10**(-5)):
		cb = []
		nodes = this.giveNodesOfSegment(segment,tol)
		cbe = np.zeros([len(nodes), 2])
		cbe[:, 0] = nodes*nv
		cbe[:, 1] = value
		cb += cbe.tolist()
		return cb

	def cbeAllBorders(this,value,tol=1*10**(-5)):
		for s in range(len(this.segments)):
			for i in range(this.nvn):
				this.cbe += this.cbFromSegment(s,value,(i+1),tol)