from FEM import NonLocalElasticityFromTensor
from FEM.Geometry import Geometry3D
from FEM.Solvers import LinealEigen
import numpy as np
import logging
from scipy.sparse.linalg import eigsh
from sendEmail import sendMailOutlook
import sys
from FEM.Utils import enmalladoEsferaFernando
import datetime


# .__class__.__name__
L = float(sys.argv[1])
l = float(sys.argv[2])
Z = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
nex = int(sys.argv[3])
omega = 6
Lr = omega*l
C11 = float(sys.argv[4])  # 223.1e6
C12 = float(sys.argv[5])  # 63.9e6
C44 = float(sys.argv[6])  # 79.6e6
rho = float(sys.argv[7])  # 2.329
C = np.array([
    [C11, C12, C12, 0.0, 0.0, 0.0],
    [C12, C11, C12, 0.0, 0.0, 0.0],
    [C12, C12, C11, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, C44, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, C44, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, C44]])


def af(rho):
    return (1/(8*np.pi*l**3))*np.exp(-rho)


coords, dicc = enmalladoEsferaFernando(2*L, nex)
geo = Geometry3D(dicc, coords, ["B1V"]*len(dicc), nvn=3, fast=True)
time = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")

log_filename = f'Sphere_{L}_{l}_{time}'

O = NonLocalElasticityFromTensor(
    geo, C, rho, l, 0.0, Lr, af, solver=LinealEigen, name=log_filename, verbose=True)
logging.info(format(sys.argv))
print("Carb")
logging.info('Creating element matrices...')
O.elementMatrices()
logging.info('Done!')
duration = O.logger.end_timer().total_seconds()
O.properties['duration'] = duration
for z in Z[::-1]:
    logging.info(f'Solving for z={z}')
    filename = f'Sphere_{L}_{l}_{z}_{time}.json'

    O.z1 = z
    O.z2 = 1-z
    O.ensembling()
    O.boundaryConditions()
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
                        secrests_path='secrets.txt', files=[f'nolocal_runner_spheres_{log_filename}.log'])
    except Exception as e:
        logging.error(e)
