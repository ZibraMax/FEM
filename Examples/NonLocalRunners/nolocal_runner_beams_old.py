from FEM import NonLocalElasticityFromTensor, NonLocalElasticity, Region1D
from FEM.Geometry import Geometry3D
from FEM.Solvers import LinealEigen
import numpy as np
import logging
from scipy.sparse.linalg import eigsh
import sys
import datetime


# .__class__.__name__
L = float(sys.argv[1])
l = float(sys.argv[2])
Z = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
nex = int(sys.argv[3])


E = float(sys.argv[4])  # 223.1e6
v = float(sys.argv[5])  # 63.9e6
rho = float(sys.argv[6])
fa = int(sys.argv[7])
h = L/10
b = L/10
kernel = 'bi-exponential'
if fa == 1:
    print('Using exponential attenuation')
    Lr = 6*l

    def af(rho):
        return (1/(8*np.pi*l**3))*np.exp(-rho)
else:
    print('Using the other kernel')
    kernel = 'Eringen-3d'
    Lr = 10*l

    def af(rho):
        xx = rho*l
        if xx < 1e-8:
            xx = 1e-8
        return 1/(4*np.pi*l*(xx)**0.5)*np.exp(-(xx/l)**0.5)


_a = L
_b = b
_c = h
nx = nex
ny = nex//10+1
nz = nex//10+1
dx = _a/nx
dy = _b/ny
dz = _c/nz

coords = []
for i in range(nx+1):
    x = i*dx
    for j in range(ny+1):
        y = j*dy
        for k in range(nz+1):
            z = k*dz
            coords += [[x, y, z]]
dicc = []


def node(i, j, k): return i*(ny+1)*(nz+1)+j*(nz+1)+k


for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            node1 = node(i, j, k)
            node2 = node(i+1, j, k)
            node3 = node(i+1, j+1, k)
            node4 = node(i, j+1, k)

            node5 = node(i, j, k+1)
            node6 = node(i+1, j, k+1)
            node7 = node(i+1, j+1, k+1)
            node8 = node(i, j+1, k+1)

            dicc += [[node1, node2, node3, node4,
                      node5, node6, node7, node8]]

geo = Geometry3D(dicc, coords, ["B1V"]*len(dicc), nvn=3, fast=True)
reg1 = Region1D(np.array([[0.0, 0.0, h/2],
                         [0.0, b, h/2]]))
reg2 = Region1D(np.array([[L, 0.0, h/2],
                         [L, b, h/2]]))
geo.addRegions([reg1, reg2])

cb = geo.cbFromRegion(0, 0.0, 1)
cb += geo.cbFromRegion(0, 0.0, 2)
cb += geo.cbFromRegion(0, 0.0, 3)

cb += geo.cbFromRegion(1, 0.0, 1)
cb += geo.cbFromRegion(1, 0.0, 2)
cb += geo.cbFromRegion(1, 0.0, 3)

time = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
log_filename = f'Beam_{L}_{l}_{time}'

O = NonLocalElasticity(
    geo, E, v, rho, l, 0.0, Lr, af, solver=LinealEigen, name=log_filename, verbose=True)
O.cbe = cb
logging.info(format(sys.argv))
logging.info('Creating element matrices...')
O.elementMatrices()
logging.info('Done!')
duration = O.logger.end_timer().total_seconds()
O.properties['duration'] = duration
for z in Z[::-1]:
    logging.info(f'Solving for z={z}')
    filename = f'Beam_{L}_{l}_{z}_{time}.json'

    O.z1 = z
    O.z2 = 1-z
    O.ensembling()
    free_dofs = O.boundaryConditions()
    K = O.K[free_dofs, :][:, free_dofs]
    M = O.M[free_dofs, :][:, free_dofs]
    K = K.tocsr()
    M = M.tocsr()
    logging.info('Solving...')
    eigv, eigvec = eigsh(
        K, 40, M, which='SM')
    idx = eigv.argsort()
    eigv = eigv[idx]
    eigvec = eigvec[:, idx]
    O.eigv = eigv
    O.eigvec = eigvec
    O.solver.solutions = []

    for i in eigvec.T:
        full_u = np.zeros(O.ngdl)
        full_u[free_dofs] = i
        O.solver.solutions.append(full_u)

    O.solver.solutions_info = [
        {'solver-type': O.solver.type, 'eigv': ei} for ei in O.eigv]
    O.solver.setSolution(0)
    logging.info('Solved!')
    O.properties['z1'] = O.z1
    O.properties['af'] = kernel

    O.exportJSON(filename)
