import matplotlib.pyplot as plt
import FEM
import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(ncols=4, nrows=3)
gss = [gs[:, :-1], gs[0, -1], gs[1, -1], gs[2, -1]]
O = FEM.importJSON('NonLocalNonHomogeneousPlate.json')
O.postProcess(gs=gss)
plt.show()
a = 0
