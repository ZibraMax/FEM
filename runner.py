from subprocess import Popen, CREATE_NEW_CONSOLE


ELES = [2, 5]
ERES = "10,15"

C11 = 223100000
C12 = 63900000
C44 = 79600000
rho = 2.329
for l in ELES:

    proc = Popen(
        f'C:/Users/ar257/Desktop/FEM/.venv/Scripts/python.exe srunner.py {ERES} {l} {C11} {C12} {C44} {rho}', creationflags=CREATE_NEW_CONSOLE)
