from FEM import PlaneStressNonLocalSparse
from FEM.Geometry import Geometry2D
from FEM.Geometry.Region import Region1D
from FEM.Utils import enmalladoFernando
import numpy as np
import logging
from sendEmail import sendMailOutlook
import sys
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import eigsh
# .__class__.__name__
L = float(sys.argv[1])
l = float(sys.argv[2])
Z = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
nex = int(sys.argv[3])
omega = 6
Lr = omega*l


E = 2.0*10**6
v = 0.2
t = L/10

dens = 2.329

# u0 = 1

l0 = 0.5/np.pi/l/l/t


def af(rho):
    return l0*np.exp(-rho)


_a = L
_b = L

nx = nex
ny = nex
coords, dicc = enmalladoFernando(_a, _b, nx, ny)

geo = Geometry2D(dicc, coords, ["C2V"]*len(dicc), nvn=2, fast=True)
regions = [Region1D([[0.0, 0.0], [0.0, L]]), Region1D([[L, 0.0], [L, L]])]
geo.addRegions(regions)
# cb = geo.cbFromRegion(0, 0.0, 1)
# cb += geo.cbFromRegion(0, 0.0, 2)
# cb += geo.cbFromRegion(1, u0, 1)

# geo.setCbe(cb)
# geo.show(draw_bc=True, label_bc=True)
# plt.show()

log_filename = f'SiPlate_disp_{L}_{l}'

O = PlaneStressNonLocalSparse(geo, E, v, t, l, 0.0, Lr=Lr,
                              af=af, rho=dens, verbose=True, name=log_filename)

logging.info('Creating element matrices...')
O.elementMatrices()
logging.info('Done!')
duration = O.logger.end_timer().total_seconds()
O.properties['duration'] = duration
for z in Z[::-1]:
    logging.info(f'Solving for z={z}')
    filename = f'SiPlate_disp_{L}_{l}_{z}.json'

    O.z1 = z
    O.z2 = 1-z
    O.F[:, :] = 0.0
    O.Q[:, :] = 0.0
    O.S[:, :] = 0.0
    O.ensembling()
    O.boundaryConditions()
    logging.info('Solving...')
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
                        secrests_path='secrets.txt', files=[f'nolocal_runner_plane_stress_{log_filename}.log'])
    except Exception as e:
        logging.error(e)
