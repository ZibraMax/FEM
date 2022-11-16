from subprocess import Popen, CREATE_NEW_CONSOLE
from threading import Thread
import sys


R = [float(i) for i in sys.argv[1].split(',')]
l = float(sys.argv[2])

i = 0
process = True


def runSim(R, l):
    global process, i
    L = R*l
    ne = max(40, min(int(L//(1.0*l)), 200))
    proc = Popen(
        f'c:/Users/david/Desktop/FEM/.venv/Scripts/python.exe nonlocal_beam_fixed_bessel.py {L} {l} {ne}', creationflags=CREATE_NEW_CONSOLE)
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
