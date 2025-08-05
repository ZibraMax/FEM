from subprocess import Popen, CREATE_NEW_CONSOLE
import pandas as pd

rho = 2.329
C11 = 166000000
C12 = 63900000
C44 = 79600000
l = 1
L = 10

print('=' * 50)
for nex in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
    print(f'Running for nex = {nex**3}')
    proc = Popen(
        f'C:/Users/ar257/Desktop/FEM/.venv/Scripts/python.exe nolocal_runner.py {L} {l} {nex} {C11} {C12} {C44} {rho}', creationflags=CREATE_NEW_CONSOLE)
