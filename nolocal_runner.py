from FEM import NonLocalElasticityFromTensor
from FEM.Geometry import Geometry3D
from FEM.Solvers import LinealEigen
import numpy as np
import logging
from scipy.sparse.linalg import eigsh
from sendEmail import sendMailOutlook
import sys


# .__class__.__name__
L = float(sys.argv[1])
l = float(sys.argv[2])
Z = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
nex = int(sys.argv[3])
omega = 6
Lr = omega*l
C11 = 223.1e6
C12 = 63.9e6
C44 = 79.6e6
C = np.array([
    [C11, C12, C12, 0.0, 0.0, 0.0],
    [C12, C11, C12, 0.0, 0.0, 0.0],
    [C12, C12, C11, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, C44, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, C44, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, C44]])
rho = 2.329


def af(rho):
    return (1/(8*np.pi*l**3))*np.exp(-rho)


_a = L
_b = L
_c = L/100

nx = nex
ny = nex
nz = max(int(nex/6), 4)

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

log_filename = f'SiPlate_{L}_{l}'

O = NonLocalElasticityFromTensor(
    geo, C, rho, l, 0.0, Lr, af, solver=LinealEigen, name=log_filename, verbose=True)

logging.info('Creating element matrices...')
O.elementMatrices()
logging.info('Done!')
duration = O.logger.end_timer().total_seconds()
O.properties['duration'] = duration
for z in Z[::-1]:
    logging.info(f'Solving for z={z}')
    filename = f'SiPlate_{L}_{l}_{z}.json'

    O.z1 = z
    O.z2 = 1-z
    O.ensembling()
    O.borderConditions()
    logging.info('Converting to csr format')
    O.K = O.K.tocsr()
    logging.info('Solving...')
    eigv, eigvec = eigsh(
        O.K, 20, O.M, which='SM')
    idx = eigv.argsort()
    eigv = eigv[idx]
    eigvec = eigvec[:, idx]
    O.eigv = eigv
    O.eigvec = eigvec
    O.solver.solutions = eigvec.T
    O.solver.solutions_info = [
        {'solver-type': O.solver.type, 'eigv': ei} for ei in O.eigv]
    O.solver.setSolution(0)
    logging.info('Solved!')

    O.exportJSON(filename)

    try:
        sendMailOutlook(mss=f"{filename} ha terminado!",
                        secrests_path='secrets.txt', files=[f'nolocal_runner_{log_filename}.log'])
    except Exception as e:
        logging.error(e)
