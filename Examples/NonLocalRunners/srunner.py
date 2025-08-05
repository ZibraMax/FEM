from subprocess import Popen, CREATE_NEW_CONSOLE
from threading import Thread
import sys


R = [float(i) for i in sys.argv[1].split(',')]
l = float(sys.argv[2])

i = 0
process = True


nex = int(sys.argv[3])
E = float(sys.argv[4])
v = float(sys.argv[5])
rho = float(sys.argv[6])
fa = int(sys.argv[7])


def runSim(R, l):
    global process, i
    L = R
    print('Running', i, 'of', len(R), 'for R =', R, 'l =',
          l, 'nex =', nex, 'E =', E, 'v =', v, 'rho =', rho)
    proc = Popen(
        f'C:/Users/ar257/Desktop/FEM/.venv/Scripts/python.exe nolocal_runner_spheres.py {L} {l} {nex} {E} {v} {rho} {fa}', creationflags=CREATE_NEW_CONSOLE)
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
