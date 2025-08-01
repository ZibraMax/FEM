from subprocess import Popen, CREATE_NEW_CONSOLE
import pandas as pd

rho = 2.329
C11 = 166000000
C12 = 63900000
C44 = 79600000
l = 1
L = 20

print('=' * 50)
for nex in [10, 12, 14, 16]:
    print(f'Running for nex = {nex}')
    proc = Popen(
        f'C:/Users/ar257/Desktop/FEM/.venv/Scripts/python.exe nolocal_runner.py {L} {l} {nex} {C11} {C12} {C44} {rho}', creationflags=CREATE_NEW_CONSOLE)
