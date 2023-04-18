import os
import tqdm
import FEM
STR_NODE = "*Node"
STR_ELEMENT = "*Element"

FILENAME = "Job.txt"
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
geo.exportJSON("AYO.json")
