import numpy as np
import matplotlib.pyplot as plt
import FEM
import pandas as pd
import matplotlib.colors as mcolors
from tqdm import tqdm

colors = list(mcolors.TABLEAU_COLORS.values())


def analitica(O, n=1000):
    l = O.l
    L = np.max(O.geometry.gdls)
    u = np.max(O.U)
    z1 = O.z1
    l0 = 1/2/l
    ebarra = u/L
    ld = -l0*(1-z1)/z1
    x = np.linspace(0.0, L, n)
    y = ebarra-ld*l/2*ebarra*(np.exp((ld*x*l-x)/l) +
                              np.exp((ld*l*L-ld*l*x-L+x)/l))
    return x, y


Z1 = 0.5
plt.figure(figsize=[15, 10])
O = FEM.importJSON('BARRA_PISANO\SiBar_pisano_50.0_0.1_0.5.json', {
                   "notCalculateNonLocal": False}, fast=True)
l = O.l
z1 = O.z1
L = np.max(O.geometry.gdls)
_X, U1, U2, U3, U = O.profile(
    [0, L/50*0.472], [L, L/50*0.472], n=1000, plot=False)
x, y = analitica(O)
plt.plot(_X, U1, '-', c=colors[0], label=f'FEM , l={l}')
plt.plot(x, y, '--', c=colors[0],
         label=f'Pisano & Fuschi (2003), l={l}')
plt.legend()
plt.grid()
plt.xlabel(r'$x$')
plt.ylabel(r'$\varepsilon_x$')
plt.tight_layout()
plt.savefig(f'L={L}.pdf')
plt.show()
