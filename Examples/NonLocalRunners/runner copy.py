from subprocess import Popen, CREATE_NEW_CONSOLE
import numpy as np

ELES = [0.1]
ERES = np.array([10])

i = 0
for l in ELES:
    for R in ERES:
        i += 1
        L = R*1.6119915*l
        l_elemento = 0.15
        ne = min(int(L//(l_elemento)), 100)
        print(ne**3, L/ne)
        proc = Popen(
            f'c:/Users/david/Desktop/FEM/.venv/Scripts/python.exe nolocal_runner.py {L} {l} {ne}', creationflags=CREATE_NEW_CONSOLE)
        if not i % 5:
            Popen.communicate(proc)
