from FEM import PlaneStressNonLocalSparse
from FEM.Geometry import Geometry2D
from FEM.Geometry.Region import Region1D
from FEM.Utils import enmalladoFernando
import numpy as np
import logging
import datetime
import sys
from scipy.sparse.linalg import eigsh
# .__class__.__name__
L = float(sys.argv[1])
l = float(sys.argv[2])
nx = int(sys.argv[3])
E = float(sys.argv[4])
v = float(sys.argv[5])
dens = float(sys.argv[6])
af = float(sys.argv[7])

Z = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
omega = 6
Lr = omega*l
h = L/10
b = L/10
t = b
l0 = 0.5/np.pi/l/l/t
att = "bi-exponential"


def af(rho):
    return l0*np.exp(-rho)


_a = L
_b = h
ny = max(nx//10, 4)
if ny % 2 != 0:
    ny += 1

coords, dicc = enmalladoFernando(_a, _b, nx, ny)

geo = Geometry2D(dicc, coords, ["C2V"]*len(dicc), nvn=2, fast=True)

cb = geo.generateBCFromCoords(0, h/2, 0, 1)
cb += geo.generateBCFromCoords(0, h/2, 0, 2)

cb += geo.generateBCFromCoords(L, h/2, 0, 1)
cb += geo.generateBCFromCoords(L, h/2, 0, 2)

geo.setCbe(cb)


time = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
log_filename = f'Beam_{L}_{l}_{nx}_{time}'
O = PlaneStressNonLocalSparse(geo, E, v, t, l, 0.0, Lr=Lr,
                              af=af, rho=dens, verbose=True, name=log_filename)
logging.info(format(sys.argv))
logging.info('Creating element matrices...')
O.elementMatrices()
logging.info('Done!')
duration = O.logger.end_timer().total_seconds()
O.properties['duration'] = duration
for z in Z[::-1]:
    logging.info(f'Solving for z={z}')
    filename = f'Beam_{L}_{l}_{nx}_{z}_{time}.json'

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
        K, 20, M, which='SM')
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
    O.properties['af'] = att

    O.exportJSON(filename)
