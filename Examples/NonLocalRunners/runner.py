from subprocess import Popen, CREATE_NEW_CONSOLE


ELES = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
ERES = "5,10,20,40,80,100"


for l in ELES:
    proc = Popen(
        f'c:/Users/david/Desktop/FEM/.venv/Scripts/python.exe srunner.py {ERES} {l}', creationflags=CREATE_NEW_CONSOLE)
