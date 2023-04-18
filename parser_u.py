import os
import numpy as np
import re

FILENAME = "U_TAPON.txt"
STR_MARK = "--------------------------------------------------------------------------------------------------"
U = []
with open(FILENAME) as f:
    line = f.readline()
    while not STR_MARK in line:
        line = f.readline()
    line = f.readline()
    while not STR_MARK in line:
        line = f.readline()
    line = f.readline()
    last = None
    while line != "\n" and line != "":
        data = re.split('\s+', line)
        n_nodo = int(data[2])
        X = float(data[4])
        Y = float(data[5])
        Z = float(data[6])
        if n_nodo != last:
            last = n_nodo
            U.append(X)
            U.append(Y)
            U.append(Z)
        line = f.readline()
    U = np.array(U)
    np.savetxt(f'{FILENAME.split(".")[0]}_PARSED.txt', U)
