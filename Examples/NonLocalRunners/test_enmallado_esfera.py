import numpy as np
import matplotlib.pyplot as plt
from FEM.Utils import enmalladoEsferaFernando
from FEM import Geometry3D

coords, con = enmalladoEsferaFernando(1, 20)

geo = Geometry3D(con, coords, ['B1V']*len(coords), 3, fast=True)
geo.exportJSON('SphereParcial.json')
