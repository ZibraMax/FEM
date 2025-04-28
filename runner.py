from subprocess import Popen, CREATE_NEW_CONSOLE


ELES = [2, 5]
ERES = "5,10,15"

C11 = 1280.0e6
C12 = 124.0e6
C44 = 578.0e6
rho = 3.515
for l in ELES:

    proc = Popen(
        f'c:/Users/david/Desktop/FEM/.venv/Scripts/python.exe srunner.py {ERES} {l} {C11} {C12} {C44} {rho}', creationflags=CREATE_NEW_CONSOLE)
