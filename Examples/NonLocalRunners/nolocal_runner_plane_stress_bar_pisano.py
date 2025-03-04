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
Z = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
nex = int(sys.argv[3])
omega = 6
Lr = omega*l


E = 2.1*10**6
v = 0.2
t = L/500

dens = 2.329

u0 = L/1000

l0 = 0.5/np.pi/l/l/t


def af(rho):
    return l0*np.exp(-rho)


_a = L
_b = L/50

nx = nex
ny = max(nex//10, 4)
coords, dicc = enmalladoFernando(_a, _b, nx, ny)

geo = Geometry2D(dicc, coords, ["C2V"]*len(dicc), nvn=2, fast=True)
regions = [Region1D([[0.0, 0.0], [0.0, _b]]), Region1D([[L, 0.0], [L, _b]])]
geo.addRegions(regions)
cb = geo.cbFromRegion(0, 0.0, 1)
cb += geo.cbFromRegion(1, u0, 1)

cb += geo.generateBCFromCoords(0.0, _b/2, 0.0, 2)

geo.setCbe(cb)
# geo.show(draw_bc=True, label_bc=True)
# plt.show()

log_filename = f'SiPlate_pisano_{L}_{l}'

O = PlaneStressNonLocalSparse(geo, E, v, t, l, 0.0, Lr=Lr,
                              af=af, verbose=True, name=log_filename)

logging.info('Creating element matrices...')
O.elementMatrices()
logging.info('Done!')
duration = O.logger.end_timer().total_seconds()
O.properties['duration'] = duration
for z in Z[::-1]:
    logging.info(f'Solving for z={z}')
    filename = f'./BARRA_PISANO/SiBar_pisano_{L}_{l}_{z}.json'

    O.z1 = z
    O.z2 = 1-z
    O.F[:, :] = 0.0
    O.Q[:, :] = 0.0
    O.S[:, :] = 0.0
    O.ensembling()
    O.boundaryConditions()
    logging.info('Solving...')
    O.K = O.K.tocsr()
    O.solver.solutions = [spsolve(O.K, O.S)]
    O.solver.solutions_info = [{'solver-type': O.solver.type}]
    O.solver.setSolution()
    logging.info('Solved!')
    O.properties['z1'] = O.z1
    O.exportJSON(filename)

    try:
        sendMailOutlook(mss=f"{filename} ha terminado!",
                        secrests_path='secrets.txt', files=[f'nolocal_runner_plane_stress_bar_pisano_{log_filename}.log'])
    except Exception as e:
        logging.error(e)
