import json
import glob
import numpy as np
from tqdm import tqdm

files = glob.glob('./PLACA_PISANO/*.json')
msg = f'L,L,t,z1,l2,v,E,d,n,filename,l/R\n'
with open('summary_placas_pisano.csv', 'w') as desy:
    desy.write(msg)
    for f in tqdm(files):
        y = json.loads(open(f).read())
        E = y['properties']['E1'][0]
        v = y['properties']['v12'][0]
        l = y['properties']['l']
        z1 = float(f.split('_')[-1][:3])
        d = y['properties']['duration']
        L = np.max(y['nodes'])
        n = len(y['dictionary'])
        t = L/10
        msg = f'{L},{L},{t},{z1},{l},{v},{E},{d},{n},{f},{l/L}\n'
        desy.write(msg)
