from subprocess import Popen, CREATE_NEW_CONSOLE
import pandas as pd

data_file = './runs.csv'
# Read the CSV file

df = pd.read_csv(data_file)
# Extract the values from the DataFrame
ELES = df['l'].unique()  # Assuming 'l' is the same for all rows

material_properties = {
    "Silicon": {
        "Density": 2.329,
        "C11": 166.0,
        "C12": 63.9,
        "C44": 79.6
    },
    "Germanium": {
        "Density": 5.323,
        "C11": 128.9,
        "C12": 48.3,
        "C44": 67.1
    },
    "Carbon": {
        "Density": 3.515,
        "C11": 1079.0,
        "C12": 124.0,
        "C44": 578.0
    }
}

for material, props in material_properties.items():
    C11 = props['C11'] * 1e6
    C12 = props['C12'] * 1e6
    C44 = props['C44'] * 1e6
    rho = props['Density']

    print('=' * 50)
    print(f"Materal {material} properties: C11 =", C11,
          "C12 =", C12, "C44 =", C44, "rho =", rho)
    for l in ELES:
        ERES = ','.join(df[df['l'] == l]['R'].unique().astype('str'))
        print("Running for l =", l, "with R =", ERES)
        proc = Popen(
            f'C:/Users/ar257/Desktop/FEM/.venv/Scripts/python.exe srunner2.py {ERES} {l} {C11} {C12} {C44} {rho}', creationflags=CREATE_NEW_CONSOLE)
