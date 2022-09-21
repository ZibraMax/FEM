import numpy as np
import matplotlib.pyplot as plt
import FEM
import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(ncols=4, nrows=3)
gss = [gs[:, :-1], gs[0, -1], gs[1, -1], gs[2, -1]]
O = FEM.importJSON('SiPlate_pisano_5.0_0.1_0.5.json')
O.postProcess(gs=gss)

_X, U1, U2, U3, U = O.profile([0, 0.019], [5, 0.019])
np.savetxt('DEFX_NUEVO_PROGRAMA.csv', np.array(
    [_X, U1]).T, delimiter=',', fmt="%s")
plt.show()
a = 0

O = FEM.importJSON('SiPlate_pisano_5.0_0.1_1.0.json')
O.postProcess(gs=gss)

_X, U1, U2, U3, U = O.profile([0, 0.019], [5, 0.019])
np.savetxt('DEFX_NUEVO_PROGRAMA_LOCAL.csv', np.array(
    [_X, U1]).T, delimiter=',', fmt="%s")
plt.show()
a = 0
