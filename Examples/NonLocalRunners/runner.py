from subprocess import Popen, CREATE_NEW_CONSOLE
import numpy as np

ELES = [0.1, 0.5, 1, 2, 5]
ERES = np.array([10, 15, 20, 50, 100])

i = 0
for l in ELES:
    for R in ERES:
        i += 1
        L = R*1.6119915*l
        ne = min(int(L//(0.1*l)), 1000)
        print(ne, l, R)
        proc = Popen(
            f'c:/Users/david/Desktop/FEM/.venv/Scripts/python.exe nolocal_runner_plane_stress_beams.py {L} {l} {ne}', creationflags=CREATE_NEW_CONSOLE)
        if not i % 5:
            Popen.communicate(proc)
