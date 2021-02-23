from .Geometry import *

class Lineal(Geometry):
	def __init__(this, lenght, n, o, nvn=1):
		this.lenght = lenght
		dictionary = []
		gdls = []
		he = this.lenght / (n)
		for i in range(0, n):
			xa = i * he
			if o==1:
				gdls += [xa]
				dictionary+=[[i,i+1]]
			else:
				gdls += [xa,xa+he/2]
				dictionary+=[[i*o,i*o+1,i*o+2]]
		gdls+=[this.lenght]
		if o==1:
			tipo = 'L1V'
		else:
			tipo = 'L2V'
		types = [tipo]*len(dictionary)
		gdls = np.array(gdls).reshape([len(gdls),1])
		Geometry.__init__(this,dictionary, gdls, types, nvn=nvn, segments=[])
