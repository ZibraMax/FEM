# %%
from scipy.linalg import eig, eigvals, eigh, eigvalsh, eigvals_banded, eig_banded
# https://scipy.github.io/devdocs/reference/generated/scipy.linalg.bandwidth.html
import scipy.linalg as lng
import logging
from scipy.sparse.linalg import spsolve, eigsh, eigs, lobpcg
import numpy as np
from FEM.Elasticity3D import ElasticityFromTensor
from FEM.Geometry.Geometry import Geometry3D
from FEM.Solvers import LinealEigen
c11 = 223.1  # MPa
c12 = 63.9  # MPa
c44 = 79.6  # MPa

C = np.array([
    [c11, c12, c12, 0, 0, 0],
    [c12, c11, c12, 0, 0, 0],
    [c12, c12, c11, 0, 0, 0],
    [0, 0, 0, c44, 0, 0],
    [0, 0, 0, 0, c44, 0],
    [0, 0, 0, 0, 0, c44]])*10**6  # Pa-3
rho = 2.329
ct = (c44/rho)**0.5

L = 20.4356
R = L/1.6119915

nx = 12
ny = 12
nz = 12

_a = L
_b = L
_c = L

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

geometria = Geometry3D(dicc, coords, ["B1V"]*len(dicc), nvn=3, fast=True)
solver = LinealEigen

# %%
O = ElasticityFromTensor(geometria, C, rho, verbose=True, solver=solver)
O.solve(path="basesolver.csv")
O.exportJSON("Cube121212.json")

# %%

# %% [markdown]
# # Sparse
# ## [eigsh](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html#scipy.sparse.linalg.eigsh):
# Simetricas reales (K,j,M,wich="SM",return_eigenvectors=True)
#
# ## [eigs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigs.html):
# No simetricas reales (K,j,M,wich="SM/SR",return_eigenvectors=False)
#
# ## [lobpcg](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lobpcg.html#scipy.sparse.linalg.lobpcg):
#
# Simetrica definida positivamente (K,x,M) #No creo quelo pueda hacer funcionar
#
# # Densas
# ## [eig](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eig.html#scipy.linalg.eig)
# No simetrica (K,M,) valores y vectores. Los saca todos
#
# ## [eigvals](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigvals.html#scipy.linalg.eigvals)
# No simetrica (K,M,) solo valores. Los saca todos
#
# ## [eigh](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html#scipy.linalg.eigh)
# Real simétrica (K,M,eigvals_only=True,subset_by_index=[0,j])
#
# ## [eigvalsh](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigvalsh.html#scipy.linalg.eigvalsh)
#
# Real simetrica (K,M,subset_by_index=[0,j]) Solo valores propios
#
# ## [eig_banded](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eig_banded.html#scipy.linalg.eig_banded)
#
# Real simetrica se necesita meter la matriz en formato de banda y es super raro. Algo asi como shape = banda + 1,M
# Segun la documentación solamente se pasan las bandas de la matriz.
#
# (omh,eigenvalues_only=True,select="i",select_range=(0,20))
#
#
# ## [eigvals_banded](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigvals_banded.html#scipy.linalg.eigvals_banded)
#
# Real simetrica se necesita meter la matriz en formato de banda y es super raro. Algo asi como shape = banda + 1,M
# Segun la documentación solamente se pasan las bandas de la matriz.
#
# (omh,select="i",select_range=(0,20))
#

# %%
geo = Geometry3D.importJSON("Cube121212.json")
O = ElasticityFromTensor(geo, C, rho, verbose=True)
# O.solve()

logging.info('Creating element matrices...')
O.elementMatrices()
logging.info('Done!')
O.ensembling()
O.borderConditions()
logging.info('Converting to csr format')
O.K = O.K.tocsr()
logging.info("Listo para lo que se venga")
omh = spsolve(O.M.tocsc(), O.K.tocsc())
logging.info("Resuelto KM")


eigv = eigsh(O.K, 20, O.M, which="SM", return_eigenvectors=False)
idx = eigv.argsort()
eigv = eigv[idx]
np.savetxt("Sparse_eigsh_kym.csv", eigv, delimiter=",", fmt="%s")
logging.info('Solved!')


# %%
geo = Geometry3D.importJSON("Cube121212.json")
O = ElasticityFromTensor(geo, C, rho, verbose=True)
# O.solve()

logging.info('Creating element matrices...')
O.elementMatrices()
logging.info('Done!')
O.ensembling()
O.borderConditions()
logging.info('Converting to csr format')
O.K = O.K.tocsr()
logging.info("Listo para lo que se venga")


eigv = eigs(O.K, 20, O.M, which="SM", return_eigenvectors=False)
idx = eigv.argsort()
eigv = eigv[idx]
np.savetxt("Sparse_eigs_kym.csv", eigv, delimiter=",", fmt="%s")
logging.info('Solved!')

# %%
geo = Geometry3D.importJSON("Cube121212.json")
O = ElasticityFromTensor(geo, C, rho, verbose=True)
# O.solve()

logging.info('Creating element matrices...')
O.elementMatrices()
logging.info('Done!')
O.ensembling()
O.borderConditions()
logging.info('Converting to csr format')
O.K = O.K.tocsr()
logging.info("Listo para lo que se venga")


eigv = eigvalsh(O.K.todense(), O.M.todense(), subset_by_index=[0, 20])
idx = eigv.argsort()
eigv = eigv[idx]
np.savetxt("Dense_eigvalsh_kym.csv", eigv, delimiter=",", fmt="%s")
logging.info('Solved!')

# %%


eigv = eigsh(omh, 20, which="SM", return_eigenvectors=False)
idx = eigv.argsort()
eigv = eigv[idx]
np.savetxt("Sparse_eigsh_mk.csv", eigv, delimiter=",", fmt="%s")
logging.info('Solved!')
logging.info("Resuelto KM")
eigv = eigs(omh, 20, which="SM", return_eigenvectors=False)
idx = eigv.argsort()
eigv = eigv[idx]
np.savetxt("Sparse_eigs_mk.csv", eigv, delimiter=",", fmt="%s")
logging.info('Solved!')
logging.info("Resuelto KM")
eigv = eigvalsh(O.K.todense(), O.M.todense(), subset_by_index=[0, 20])
idx = eigv.argsort()
eigv = eigv[idx]
np.savetxt("Dense_eigvalsh_mk.csv", eigv, delimiter=",", fmt="%s")
logging.info('Solved!')
