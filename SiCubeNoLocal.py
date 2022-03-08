import numpy as np
from FEM.Elasticity3D import NonLocalElasticityFromTensor
from FEM.Geometry.Geometry import Geometry3D
from FEM.Solvers import LinealEigen
from FEM import logging

for l in [0.2, 0.4, 0.535, 0.6, 0.8, 1]:
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

    nx = 17
    ny = 17
    nz = 17

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
    z1 = 0.5
    Lr = 9*l
    print(l, Lr)

    def af(rho):
        return (1/(8*np.pi*l**3))*np.exp(-rho)  # No referencia, sacada a mano

    O = NonLocalElasticityFromTensor(
        geometria, C, rho, l, z1, Lr, af, verbose=True, solver=solver)
    O.solve(path=f"SiCube_{l}.csv")
    O.exportJSON(f"SiCube_{l}.json")
    logging.info(f"l={l} terminado")
