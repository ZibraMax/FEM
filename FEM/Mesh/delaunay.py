import triangle as tr

import numpy as np
import matplotlib.pyplot as plt

from .Geometry import *

class Delaunay1V(Geometry):
	def __init__(this, vertices, params, nvn=1):
		seg = []
		for i in range(len(vertices)-1):
			seg.append([i,i+1])
		seg.append([i+1,0])
		original = dict(vertices=np.array(vertices),segments=np.array(seg))
		triangular = tr.triangulate(original,params)
		dictionary = triangular['triangles'].tolist()
		tipos = np.zeros([len(dictionary)]).astype(str)
		if 'o2' in params:
			tipos[:] = 'T2V'
		else:
			tipos[:] = 'T1V'
		gdls = triangular['vertices'].tolist()
		if tipos[0] == 'T2V':
			for dicc in dictionary:
					a1 = dicc[5]
					a2 = dicc[3]
					a3 = dicc[4]
					dicc[3] = a1
					dicc[4] = a2
					dicc[5] = a3
		Geometry.__init__(this,dictionary,gdls,tipos,nvn=nvn,segments=seg)
		this.mask = vertices

def _strdelaunay(constrained=True,delaunay=True,a=None,q=None,o=1):
	p = ''
	if o==2:
		o = '-o2'
	else:
		o = ''
	if constrained:
		p = 'p'
	if a == None:
		a = ''
	else:
		a = 'a'+format(a)
	D = ''
	if delaunay:
		D = 'D'
	if q == None:
		q=''
	else:
		if type(q) == int:
			if q > 35:
				raise "No se puede crear una triangulacion con angulos menores a 35 grados"
		q = 'q'+format(q)
	return p+a+D+q+'i'+o