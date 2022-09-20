from FEM import PlaneStressNonLocalSparse
from FEM.Geometry import Geometry2D
from FEM.Utils import enmalladoFernando
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


E = 2.1*10**6
v = 0.2
dens = 2.329


_a = L
_b = L/10
t = L/10

l0 = 0.5/np.pi/l/l/t


def af(rho):
    return np.exp(-rho)


nx = nex
ny = min(nex//10, 4)
coords, dicc = enmalladoFernando(_a, _b, nx, ny)

geo = Geometry2D(dicc, coords, ["C2V"]*len(dicc), nvn=2, fast=True)

log_filename = f'SiBeam_{L}_{l}'

O = PlaneStressNonLocalSparse(geo, E, v, t, l, 0.0, Lr=Lr,
                              af=af, rho=dens, verbose=True, solver=LinealEigen, name=log_filename)

logging.info('Creating element matrices...')
O.elementMatrices()
logging.info('Done!')
duration = O.logger.end_timer().total_seconds()
O.properties['duration'] = duration
for z in Z[::-1]:
    logging.info(f'Solving for z={z}')
    filename = f'SiBeam_{L}_{l}_{z}.json'

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
    O.properties['z1'] = O.z1
    O.exportJSON(filename)

    try:
        sendMailOutlook(mss=f"{filename} ha terminado!",
                        secrests_path='secrets.txt', files=[f'nolocal_runner_plane_stress_beams_{log_filename}.log'])
    except Exception as e:
        logging.error(e)
