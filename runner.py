from subprocess import Popen, CREATE_NEW_CONSOLE
import numpy as np

ELES = [0.1, 0.5, 1, 2, 5]
ERES = np.array([100])

i = 0
for l in ELES:
    for R in ERES:
        i += 1
        L = R*1.6119915*l
        ne = min(int(L//(1.5*l)), 30)
        print(ne)
        proc = Popen(
            f'c:/Users/david/Desktop/FEM/.venv/Scripts/python.exe nolocal_runner.py {L} {l} {ne}', creationflags=CREATE_NEW_CONSOLE)
        if not i % 1:
            Popen.communicate(proc)
