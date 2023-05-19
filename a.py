import numpy as np
import FEM
from tqdm import tqdm

F1 = FEM.importJSON("Examples/Mesh_tests/AYO_TAPON.json", {}, fast=True)
F2 = FEM.importJSON("Examples/Mesh_tests/AYO_NO_TAPON.json", {}, fast=True)


def error(LDUM, HDUM, nvn=3):
    errors = []
    for i in tqdm(range(len(HDUM.geometry.gdls))):
        coord = HDUM.geometry.gdls[i]
        potential = LDUM.geometry.KDTree.query(coord, 40, workers=-1)[1]
        flag = True
        errores = []
        for eidx in potential:
            e = LDUM.elements[eidx]
            z = e.inverseMapping(coord.reshape([3, 1]))
            zz = np.sum(z)
            errores.append(f"{z.flatten()}, {zz}")
            tol = 1e-6
            if zz <= (1+tol) and np.all(z.flatten() >= -tol):
                uLDUMx, uLDUMy, uLDUMz = e.giveSolutionPoint(z)[1]
                uHDUMx = HDUM.U[nvn*i]
                uHDUMy = HDUM.U[nvn*i + 1]
                uHDUMz = HDUM.U[nvn*i + 2]
                errors.append(abs(uLDUMx-uHDUMx))
                errors.append(abs(uLDUMy-uHDUMy))
                errors.append(abs(uLDUMz-uHDUMz))
                flag = False
                break
        if flag:
            print(f"HEY coord {i}, {coord} NOT WORKED! {errores}")
            break
    return np.array(errors)


U = error(F1, F2)
np.savetxt("errores.txt", U)
print()
