import numpy as np
import matplotlib.pyplot as plt
import FEM
from FEM import Mesh

G = 1
phi =1

dicc = np.loadtxt('web_elements.txt',dtype=int)-1
gdls = np.loadtxt('web_nodes.txt')
types = ['C1V']*len(dicc)
segments = [[0,9],[9,19],[19,29],[29,39],[39,49],[49,59],[59,69],[69,79],[79,89],[89,99],[99,90],[90,80],[80,70],[70,60],[60,50],[50,40],[40,30],[30,20],[20,10],[10,0]]
geometria = Mesh.Geometry(dicc, gdls, types, nvn=1, segments=segments)
geometria.saveMesh('Web_test')
geometria.show()
plt.show()
O = FEM.Torsion2D(geometria,G,phi)
O.solve()
plt.show()
