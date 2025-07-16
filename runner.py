from subprocess import Popen, CREATE_NEW_CONSOLE
import pandas as pd

data_file = './runs.csv'
df = pd.read_csv(data_file)
ELES = df['l'].unique()


print('=' * 50)
for l in ELES:
    ERES = ','.join(df[df['l'] == l]['R'].unique().astype('str'))
    print("Running for l =", l, "with R =", ERES)
    nex = 30
    E = 30e6
    v = 0.3
    rho = 1
    fa = 1
    proc = Popen(
        f'C:/Users/ar257/Desktop/FEM/.venv/Scripts/python.exe srunner2.py {ERES} {l} {nex} {E} {v} {rho} {fa}', creationflags=CREATE_NEW_CONSOLE)
