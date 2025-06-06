from subprocess import Popen, CREATE_NEW_CONSOLE


ELES = [1, 2, 4, 6]
ERESS = {"1": "2.5,1.25",
         "2": "5,2.5",
         "4": "10,5",
         "6": "15,7.5", }

C11 = 182500000
C12 = 48300000
C44 = 67100000
rho = 5.323
for l in ELES:
    ERES = ERESS[str(l)]
    proc = Popen(
        f'C:/Users/ar257/Desktop/FEM/.venv/Scripts/python.exe srunner2.py {ERES} {l} {C11} {C12} {C44} {rho}', creationflags=CREATE_NEW_CONSOLE)
