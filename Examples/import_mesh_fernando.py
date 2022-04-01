from FEM.Geometry import Geometry3D

# SPHERE MESH

gdls = []
els = []
with open('input.txt') as f:
    ngdl, nele = [int(i) for i in f.readline().split('\t')][:2]
    for i in range(ngdl):
        x, y, z = [float(i) for i in f.readline().split('\t')]
        gdls.append([x, y, z])
    for i in range(nele):
        ele = [int(i)-1 for i in f.readline().split('\t')][:8]
        els.append(ele)

types = ['B1V']*len(els)
geo = Geometry3D(els, gdls, types, nvn=3, fast=True)
geo.exportJSON('SPHERE_FERNANDO_161616.json')

# # PIRAMID MESH

# gdls = []
# els = []
# with open('mesh_piramid.txt') as f:
#     ngdl, nele = [int(i) for i in f.readline().split('\t')][:2]
#     for i in range(ngdl):
#         x, y, z = [float(i) for i in f.readline().split('\t')]
#         gdls.append([x, y, z])
#     for i in range(nele):
#         ele = [int(i)-1 for i in f.readline().split('\t')][:8]
#         els.append(ele)

# types = ['B1V']*len(els)
# geo = Geometry3D(els, gdls, types, nvn=3, fast=True)
# geo.exportJSON('PIRAMID_FERNANDO.json')
