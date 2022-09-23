import numpy as np
import matplotlib.pyplot as plt
import FEM
import pandas as pd
import matplotlib.colors as mcolors
from tqdm import tqdm

colors = list(mcolors.TABLEAU_COLORS.values())

df = pd.read_csv('summary_placas_pisano.csv')
df.head()

LS = df['l/R'].unique()

Z1 = 0.5
plt.figure(figsize=[10, 8])
for L in tqdm(LS):
    df_chiquito = df[df["z1"] == Z1]
    df_chiquito = df_chiquito[df_chiquito["l/R"] == L]
    for i in range(len(df_chiquito)):
        O = FEM.importJSON(
            df_chiquito['filename'].values[i], {"notCalculateNonLocal": False}, fast=True)
        l = O.l
        z1 = O.z1
        L = np.max(O.geometry.gdls)
        print(L, l)
        _X, U1, U2, U3, U = O.profile(
            [0, 0.019], [L, 0.019], n=1000, plot=False)
        plt.plot(np.array(_X)/L, U1, '-', c=colors[i], label=f'FEM , l={l}')
    plt.legend()
    plt.grid()
    plt.xlabel(r'$x/L$')
    plt.ylabel(r'$\varepsilon_x$')
    plt.tight_layout()
    plt.savefig(f'l_R={l/L:.3f}.pdf')
    plt.show()
