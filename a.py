import numpy as np
import FEM
F1 = FEM.importJSON("Examples/Mesh_tests/AYO_TAPON.json",
                    {"verbose": True}, fast=True)
F2 = FEM.importJSON("Examples/Mesh_tests/AYO_NO_TAPON.json",
                    {"verbose": True}, fast=True)


def error(LDUM, HDUM, nvn=3):
    errors = []
    for i in range(len(HDUM.geometry.gdls)):
        coord = HDUM.geometry.gdls[i]
        for e in LDUM.elements:
            if e.isInside(coord.reshape([1, 3]))[0]:
                z = e.inverseMapping(coord.reshape([3, 1]))
                uLDUMx, uLDUMy, uLDUMz = e.giveSolutionPoint(z)[1]
                uHDUMx = HDUM.U[nvn*i]
                uHDUMy = HDUM.U[nvn*i + 1]
                uHDUMz = HDUM.U[nvn*i + 2]
                errors.append(abs(uLDUMx-uHDUMx))
                errors.append(abs(uLDUMy-uHDUMy))
                errors.append(abs(uLDUMz-uHDUMz))
                break
    return np.array(errors)


U = error(F1, F2)
np.savetxt("errores.txt", U)
print()
