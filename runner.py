from subprocess import Popen, CREATE_NEW_CONSOLE
import pandas as pd

data_file = './runs.csv'
df = pd.read_csv(data_file)

print('=' * 50)
for row in df.iloc:
    R = row['R']
    l = row['l']
    nex = int(row['nex'])
    E = row['E']
    v = row['v']
    rho = row['rho']
    af = row['af']
    k = row['k']

    print("Running with parameters:")
    print(
        f"R: {R}, l: {l}, nex: {nex}, E: {E}, v: {v}, rho: {rho}, af: {af}, k: {k}")
    proc = Popen(
        f'C:/Users/ar257/Desktop/FEM/.venv/Scripts/python.exe nolocal_runner_beams.py {R} {l} {nex} {E} {v} {rho} {af} {k}', creationflags=CREATE_NEW_CONSOLE)
