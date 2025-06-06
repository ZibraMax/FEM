import numpy as np
from subprocess import Popen, CREATE_NEW_CONSOLE
from threading import Thread
import sys


R = [float(i) for i in sys.argv[1].split(',')]
l = float(sys.argv[2])

i = 0
process = True


C11 = float(sys.argv[3])  # 223.1e6
C12 = float(sys.argv[4])  # 63.9e6
C44 = float(sys.argv[5])  # 79.6e6
rho = float(sys.argv[6])  # 2.329


def runSim(R, l):
    global process, i
    L = R*1.6119915
    ne = max(7, min(int(L//(l)), 200))
    print(L, l, ne)
    proc = Popen(
        f'C:/Users/ar257/Desktop/FEM/.venv/Scripts/python.exe nolocal_runner.py {L} {l} {ne} {C11} {C12} {C44} {rho}', creationflags=CREATE_NEW_CONSOLE)
    Popen.communicate(proc)
    i = i + 1
    process = False


def check():
    global process, i, l, R
    while process:
        print('Running', i)
        # here you can wait as much as you want without freezing the program
        runSim(R[i], l)
    else:
        print('A')
        if i < len(R):
            process = True
            check()
        print('Finished')


t = Thread(target=check)
t.deamon = True
t.start()
