import numpy as np
import os
import tqdm
import FEM
STR_NODE = "*Node"
STR_ELEMENT = "*Element"

FILENAME = "NO_TAPON.txt"
nodes = []
elements = []

with open(FILENAME, encoding="utf-8") as f:
    line = f.readline()
    while not STR_NODE in line:
        line = f.readline()
    pb = tqdm.tqdm()
    line = f.readline()
    while not STR_ELEMENT in line:
        nodes.append([float(i) for i in line.split(",")[1:]])
        line = f.readline()
        pb.update()
    pb = tqdm.tqdm()
    line = f.readline()
    while not "*" in line and line != "":
        cone = [int(i)-1 for i in line.split(",")[1:]]
        elements.append(cone)
        line = f.readline()
        pb.update()
nel = len(elements)
geo = FEM.Geometry3D(elements, nodes, ["TE2V"]*nel, 3, fast=True)
O = FEM.Elasticity(geo, 1, 0.3, 0.1)
U = np.loadtxt(f"U_{FILENAME.split('.')[0]}_PARSED.txt")
O.solver.solutions = [U]
O.solver.solutions_info = [{'solver-type': "ABAQUS!"}]
O.solver.setSolution()
O.exportJSON(f"AYO_{FILENAME.split('.')[0]}.json")
