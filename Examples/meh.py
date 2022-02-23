import json
dicc = []
coords = []
with open('mesh.csv', 'r') as f:
    nnodes, nels, ngdl = [int(i) for i in f.readline().split(',')]
    for _ in range(nels):
        dicc += [[int(i)-1 for i in f.readline().split(',')][1:]]
    for _ in range(ngdl):
        coords += [[float(i) for i in f.readline().split(',')][1:]]
types = ['B1V']*nels

x = {
    "nodes": coords,
    "dictionary": dicc,
    "types": types,
    "nvn": 3,
    "ngdl": ngdl,
    "disp_field": [[0.0]*3*ngdl]
}
y = json.dumps(x)
with open('../FEM-C--/docs/resources/piramid.json', "w") as f:
    f.write(y)

a = 0
