from subprocess import Popen, CREATE_NEW_CONSOLE


ELES = [1]
ERES = "10"


for l in ELES:
    proc = Popen(
        f'c:/Users/david/Desktop/FEM/.venv/Scripts/python.exe srunner.py {ERES} {l}', creationflags=CREATE_NEW_CONSOLE)
